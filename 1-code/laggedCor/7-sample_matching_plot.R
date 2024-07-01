calculate_match_index <-
  function(x,
           y,
           time1,
           time2,
           time_tol = 5,
           step = 0.5,
           threads = 10) {
    time_window1 <-
      seq(from = step / 2, to = time_tol, by = step)
    
    time_window2 <-
      -rev(seq(from = step / 2, to = time_tol, by = step))
    
    time_window <- sort(c(time_window2, time_window1))
    
    shift_time <-
      paste("(",
            paste(round(time_window[-length(time_window)] * 60, 2),
                  round(time_window[-1] * 60, 2), sep = ','),
            "]", sep = "")
    
    temp_fun <-
      function(temp_idx,
               time_window,
               x,
               y,
               time1,
               time2) {
        idx =
          time1 %>%
          purrr::map(function(x) {
            diff_time =
              difftime(x, time2, units = "hours")
            which(diff_time > time_window[temp_idx] &
                    diff_time <= time_window[temp_idx + 1])
          })
      }
    
    bpparam <-
      BiocParallel::MulticoreParam(workers = threads,
                                   progressbar = TRUE)
    
    idx <-
      BiocParallel::bplapply(
        X = 1:(length(time_window) - 1),
        FUN = temp_fun,
        time_window = time_window,
        x = x,
        y = y,
        time1 = time1,
        time2 = time2,
        BPPARAM = bpparam
      )
    
    shift_time_mean <-
      lapply(
        1:(length(time_window) - 1),
        FUN = function(idx) {
          mean(time_window[c(idx, idx + 1)])
        }
      ) %>%
      unlist() %>%
      `*`(60)
    
    
    all_idx <-
      list(
        idx = idx,
        time_window = time_window,
        shift_time = shift_time,
        shift_time_mean = shift_time_mean
      )
    return(all_idx)
  }


match_index_plot1 <-
  function(object,
           ...) {
    UseMethod("match_index_plot1")
  }

match_index_plot1.list <-
  function(object,
           ....) {
    plot <-
      object$idx %>%
      lapply(function(x) {
        length(unique(unlist(x)))
      }) %>%
      unlist() %>%
      data.frame(number = .) %>%
      dplyr::mutate(index = object$shift_time_mean) %>%
      ggplot(aes(index, number)) +
      geom_segment(aes(
        x = index,
        y = 160,
        xend = index,
        yend = number
      )) +
      geom_point(size = 3) +
      base_theme +
      geom_text(aes(x = index, y = number, label = number),
                nudge_x = 12,
                nudge_y = 4) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "Shift time (min)", y = "Matched sample number")
    
    return(plot)
    
  }

match_index_plot1.lagged_cor_result <-
  function(object,
           ....) {
    object <- object@idx
    plot <-
      object$idx %>%
      lapply(function(x) {
        length(unique(unlist(x)))
      }) %>%
      unlist() %>%
      data.frame(number = .) %>%
      dplyr::mutate(index = object$shift_time_mean) %>%
      ggplot(aes(index, number)) +
      geom_segment(aes(
        x = index,
        y = 160,
        xend = index,
        yend = number
      )) +
      geom_point(size = 3) +
      base_theme +
      geom_text(aes(x = index, y = number, label = number),
                nudge_x = 12,
                nudge_y = 4) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "Shift time (min)", y = "Matched sample number")
    
    return(plot)
    
  }


match_index_plot2 <-
  function(object,
           ...) {
    UseMethod("match_index_plot2")
  }

match_index_plot2.list <-
  function(object,
           ....) {
    temp_data <-
      purrr::map2(
        .x = object$idx,
        .y = object$shift_time_mean,
        .f = function(idx,
                      t) {
          purrr::map2(
            time1,
            idx,
            .f = function(x, y) {
              difftime(time1 = x,
                       time2 = time1[y],
                       units = "hour")
            }
          ) %>%
            unlist() %>%
            data.frame(time = ., shift_time = t) %>%
            dplyr::mutate(time = time * 60)
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    plot <-
      temp_data %>%
      dplyr::mutate(shift_time = as.character(shift_time)) %>%
      dplyr::mutate(shift_time = factor(shift_time,
                                        levels = as.character(unique(
                                          temp_data$shift_time
                                        )))) %>%
      ggplot(aes(x = shift_time, time)) +
      geom_half_boxplot(
        aes(x = shift_time, time),
        fill = "transparent",
        show.legend = FALSE,
        side = "l"
      ) +
      geom_half_violin(
        aes(x = shift_time, time, color = shift_time),
        show.legend = FALSE,
        side = "r"
      ) +
      geom_half_point(
        aes(x = shift_time, time, color = shift_time),
        show.legend = FALSE,
        side = "l",
        alpha = 0.7
      ) +
      ggsci::scale_color_lancet() +
      base_theme +
      labs(y = "Different time (min)", x = "Shift time (min)")
    
    return(plot)
    
  }

match_index_plot1.lagged_cor_result <-
  function(object,
           ....) {
    object <- object@idx
    
    temp_data <-
      purrr::map2(
        .x = object$idx,
        .y = object$shift_time_mean,
        .f = function(idx,
                      t) {
          purrr::map2(
            time1,
            idx,
            .f = function(x, y) {
              difftime(time1 = x,
                       time2 = time1[y],
                       units = "hour")
            }
          ) %>%
            unlist() %>%
            data.frame(time = ., shift_time = t) %>%
            dplyr::mutate(time = time * 60)
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    plot <-
      temp_data %>%
      dplyr::mutate(shift_time = as.character(shift_time)) %>%
      dplyr::mutate(shift_time = factor(shift_time,
                                        levels = as.character(unique(
                                          temp_data$shift_time
                                        )))) %>%
      ggplot(aes(x = shift_time, time)) +
      gghalves::geom_half_boxplot(
        aes(x = shift_time, time),
        fill = "transparent",
        show.legend = FALSE,
        side = "l"
      ) +
      gghalves::geom_half_violin(
        aes(x = shift_time, time, color = shift_time),
        show.legend = FALSE,
        side = "r"
      ) +
      gghalves::geom_half_point(
        aes(x = shift_time, time, color = shift_time),
        show.legend = FALSE,
        side = "l",
        alpha = 0.7
      ) +
      ggsci::scale_color_lancet() +
      base_theme +
      labs(y = "Different time (min)", x = "Shift time (min)")
    
    return(plot)
    
    
  }




sample_matching_plot <-
  function(object,
           index = 1,
           only_remain_matched = TRUE,
           day_night_df = NULL,
           add_text = FALSE) {
    match_idx = object@idx[[index]]
    time1 = object@time1
    time2 = object@time2
    
    idx1 =
      lapply(match_idx, length) %>%
      unlist() %>%
      `!=`(0) %>%
      which()
    
    if (length(idx1) == 0) {
      return(NULL)
    }
    
    idx2 = match_idx[idx1]
    
    segment_data =
      purrr::map2(.x = idx1, .y = idx2, function(x, y) {
        data.frame(
          x = time1[x],
          xend = time2[y],
          y = "time1",
          yend = "time2"
        )
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    if (only_remain_matched) {
      temp_data =
        data.frame(x = c(time1[unique(idx1)], time2[unique(unlist(idx2))]),
                   y = c(rep("time1", length(time1[unique(idx1)])),
                         rep("time2", length(time2[unique(unlist(idx2))]))))
    } else{
      temp_data =
        data.frame(x = c(time1, time2),
                   y = c(rep("time1", length(time1)), rep("time2", length(time2))))
    }
    
    if (!is.null(day_night_df)) {
      plot =
        ggplot() +
        geom_rect(
          mapping = aes(
            xmin = start,
            xmax = end,
            ymin = -Inf,
            ymax = Inf
          ),
          fill = "lightyellow",
          data = day_night_df,
          show.legend = FALSE
        ) +
        geom_segment(aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ),
        data = segment_data) +
        geom_point(aes(x = x, y = y,
                       color = y),
                   data = temp_data,
                   show.legend = FALSE) +
        scale_x_datetime(
          breaks = scales::date_breaks("12 hour"),
          date_labels = "%a %H:%M",
          timezone = "America/Los_Angeles"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        panel.grid = element_blank()) +
        labs(x = "", y = "") +
        ggsci::scale_color_jama()
    } else{
      plot =
        ggplot() +
        geom_segment(aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ),
        data = segment_data) +
        geom_point(aes(x = x, y = y,
                       color = y),
                   data = temp_data,
                   show.legend = FALSE) +
        scale_x_datetime(
          breaks = scales::date_breaks("12 hour"),
          date_labels = "%a %H:%M",
          timezone = "America/Los_Angeles"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        panel.grid = element_blank()) +
        labs(x = "", y = "") +
        ggsci::scale_color_jama()
    }
    
    shift_time =
      object$shift_time[index]
    
    plot =
      plot +
      ggtitle(
        label = paste("Shift time(min):", shift_time, sep = ""),
        subtitle = paste(length(idx1), "matched samples")
      )
    
    if (add_text) {
      plot =
        plot +
        ggrepel::geom_text_repel(
          mapping = aes(
            x = x,
            y = y,
            label = paste(lubridate::hour(x), lubridate::minute(x), sep = ":")
          ),
          data = temp_data,
          size = 2
        )
    }
    
    plot
    
  }
