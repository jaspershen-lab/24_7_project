no_function()

library(tidyverse)
setwd(masstools::get_project_wd())
rm(list = ls())
source("1_code/tools.R")

###internal omics data
{
  ###cortisol
  load("3_data_analysis/data_preparation/cortisol/sample_info")
  load("3_data_analysis/data_preparation/cortisol/variable_info")
  cortisol_sample_info = sample_info
  cortisol_variable_info = variable_info
  
  ###cytokine
  load("3_data_analysis/data_preparation/cytokine/sample_info")
  load("3_data_analysis/data_preparation/cytokine/variable_info")
  cytokine_sample_info = sample_info
  cytokine_variable_info = variable_info
  
  ###lipidomics
  load("3_data_analysis/data_preparation/lipidomics/sample_info")
  load("3_data_analysis/data_preparation/lipidomics/variable_info")
  lipidomics_sample_info = sample_info %>%
    dplyr::rename(accurate_time = DT)
  lipidomics_variable_info = variable_info
  
  ###metabolic_panel
  load("3_data_analysis/data_preparation/metabolic_panel/sample_info")
  load("3_data_analysis/data_preparation/metabolic_panel/variable_info")
  metabolic_panel_sample_info = sample_info
  metabolic_panel_variable_info = variable_info
  
  ###metabolomics
  load("3_data_analysis/data_preparation/metabolomics/sample_info")
  load("3_data_analysis/data_preparation/metabolomics/variable_info")
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  
  ##remove QC and blank samples
  metabolomics_sample_info =
    metabolomics_sample_info %>%
    dplyr::filter(!is.na(accurate_time))
  
  ###proteomics
  load("3_data_analysis/data_preparation/proteomics/sample_info")
  load("3_data_analysis/data_preparation/proteomics/variable_info")
  proteomics_sample_info = sample_info
  proteomics_variable_info = variable_info
  
  ###total protein
  load("3_data_analysis/data_preparation/total_protein/sample_info")
  load("3_data_analysis/data_preparation/total_protein/variable_info")
  total_protein_sample_info = sample_info
  total_protein_variable_info = variable_info
}

###set work directory
setwd("3_data_analysis/summary_info")

library(scales)

######internal omics data
temp_data <-
  rbind(
    data.frame(time = total_protein_sample_info[, c("accurate_time")],
               class = "total_protein"),
    data.frame(time = cortisol_sample_info[, c("accurate_time")],
               class = "cortisol"),
    data.frame(time = cytokine_sample_info[, c("accurate_time")],
               class = "cytokine"),
    data.frame(time = lipidomics_sample_info[, c("accurate_time")],
               class = "lipidomics"),
    data.frame(time = metabolic_panel_sample_info[, c("accurate_time")],
               class = "metabolic_panel"),
    data.frame(time = metabolomics_sample_info[, c("accurate_time")],
               class = "metabolomics"),
    data.frame(time = proteomics_sample_info[, c("accurate_time")],
               class = "proteomics")
  ) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(time)) %>%
  dplyr::mutate(class = factor(class, levels = unique(class)))

sun_rise_time = "6:00:00"
sun_set_time = "18:00:00"

sun_rise <-
  lubridate::ymd_hms(paste(unique(lubridate::date(temp_data$time)), c(sun_rise_time)),
                     tz = lubridate::tz(temp_data$time))

sun_set =
  lubridate::ymd_hms(paste(unique(lubridate::date(temp_data$time)), c(sun_set_time)),
                     tz = lubridate::tz(temp_data$time))

day_df <-
  data.frame(start = sun_rise,
             end = sun_set,
             day = lubridate::date(sun_rise)) %>%
  dplyr::mutate(week = format(day, "%a")) %>%
  dplyr::mutate(week = paste(week,
                             lubridate::month(day),
                             lubridate::day(day),
                             sep = "-")) %>%
  dplyr::mutate(week = factor(week, unique(week)))

min(temp_data$time)
min_time = lubridate::as_datetime("2019-04-29 03:00:00", tz = "America/Los_Angeles")
max(temp_data$time)
max_time = lubridate::as_datetime("2019-05-07 07:00:00", tz = "America/Los_Angeles")

head(day_df)
tail(day_df)

day_df$end[nrow(day_df)] <- max_time + 60 * 60

night_df <-
  data.frame(start = day_df$end[-nrow(day_df)],
             end = day_df$start[-1])

night_df <-
  rbind(data.frame(start = min_time - 60 * 60,
                   end = day_df$start[1]),
        night_df)

plot1 <-
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_df,
    show.legend = FALSE
  ) +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = alpha("grey", 0.2),
    data = night_df,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = time, y = class, fill = class),
    color = "black",
    show.legend = FALSE,
    shape = 21,
    size = 1,
    data = temp_data
  ) +
  scale_fill_manual(values = class_color) +
  scale_x_datetime(
    expand = expansion(mult = c(0, 0)),
    breaks = scales::date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    timezone = "America/Los_Angeles"
  ) +
  labs(x = "Time", y = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 10
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot2 =
  temp_data %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(class = factor(class, levels = class)) %>%
  ggplot(aes(x = n, y = class)) +
  geom_bar(stat = "identity", aes(fill = class), show.legend = FALSE) +
  labs(x = "Sample number", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = class_color) +
  geom_text(aes(x = n + 2, y = class, label = n)) +
  geom_text(aes(x = n + 2, y = class, label = n)) +
  base_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot2

library(patchwork)

plot =
  plot1 + plot2 + patchwork::plot_layout(ncol = 2, widths = c(10, 1))
plot
ggsave(plot,
       file = "internal_omics.pdf",
       width = 14,
       height = 2)

######sample collection overiview
temp_data =
  rbind(
    data.frame(time = total_protein_sample_info[, c("accurate_time")],
               class = "total_protein"),
    data.frame(time = cortisol_sample_info[, c("accurate_time")],
               class = "cortisol"),
    data.frame(time = cytokine_sample_info[, c("accurate_time")],
               class = "cytokine"),
    data.frame(time = lipidomics_sample_info[, c("accurate_time")],
               class = "lipidomics"),
    data.frame(time = metabolic_panel_sample_info[, c("accurate_time")],
               class = "metabolic_panel"),
    data.frame(time = metabolomics_sample_info[, c("accurate_time")],
               class = "metabolomics"),
    data.frame(time = proteomics_sample_info[, c("accurate_time")],
               class = "proteomics")
  )

library(lubridate)
library(hms)
accurate_time = unique(temp_data$time)

time =
  ymd_hms(accurate_time) %>%
  hms::as_hms() %>%
  as.numeric()
# as.POSIXct()

time = time / 60 / 60

temp_data =
  data.frame(accurate_time, time) %>%
  dplyr::mutate(day = lubridate::date(accurate_time),
                hour = lubridate::hour(accurate_time)) %>%
  dplyr::arrange(day, time) %>%
  dplyr::mutate(image = "blood.png")

library(ggimage)

plot <-
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_df,
    show.legend = FALSE
  ) +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = alpha("grey", 0.2),
    data = night_df,
    show.legend = FALSE
  ) +
  geom_image(
    aes(x = accurate_time,
        y = time,
        image = image),
    asp = 2,
    size = 0.05,
    data = temp_data
  ) +
  scale_x_datetime(
    expand = expansion(mult = c(0, 0)),
    breaks = scales::date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    timezone = "America/Los_Angeles"
  ) +
  labs(x = "Time", y = "") +
  scale_y_continuous(breaks = seq(0, 24, 2),
                     labels = seq(0, 24, 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 10
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot
plot <-
  plot + plot2 + patchwork::plot_layout(ncol = 2, widths = c(10, 1))

ggsave(plot,
       filename = "sample_collection.pdf",
       width = 14,
       height = 3)

#######sampling frequency
library(plyr)

day_class_info =
  day_df %>%
  dplyr::mutate(class = 1:9)

class =
  temp_data$accurate_time %>%
  purrr::map(function(x) {
    class = which(x > day_class_info$start &
                    x < day_class_info$end)
    if (length(class) == 0) {
      class = 0
    }
    class
  }) %>%
  unlist()

temp_data$class = class

library(gghalves)

temp_data2 =
  temp_data %>%
  dplyr::filter(class != 0) %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::arrange(time) %>%
      dplyr::pull(time) %>%
      diff()
  }) %>%
  unlist() %>%
  data.frame(class = "yes", time = .)

plot =
  temp_data2 %>%
  ggplot(aes(x = class, y = time)) +
  labs(x = "", y = "Sampling interval (hour)") +
  geom_half_violin() +
  geom_half_boxplot(
    fill = "transparent",
    color = ggsci::pal_aaas()(n = 10)[2],
    outlier.shape = NA
  ) +
  geom_dotplot(binaxis = "y",
               method = "histodot",
               stackdir = "up") +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot

ggsave(plot, filename = "sampling_frequency.pdf", width = 3, height = 7)

quantile(temp_data2$time)

temp_data =
  data.frame(
    number = c(
      nrow(total_protein_variable_info),
      nrow(cortisol_variable_info),
      nrow(cytokine_variable_info),
      nrow(lipidomics_variable_info),
      nrow(metabolic_panel_variable_info),
      nrow(metabolomics_variable_info),
      nrow(proteomics_variable_info)
    ),
    class = c(
      "total_protein",
      "cortisol",
      "cytokine",
      "lipidomics",
      "metabolic_panel",
      "metabolomics",
      "proteomics"
    )
  ) %>%
  dplyr::mutate(class = factor(class, levels = class))

plot =
  temp_data %>%
  ggplot(aes(x = number, y = class)) +
  geom_segment(aes(
    x = 0,
    y = class,
    xend = number,
    yend = class,
    color = class
  ),
  show.legend = FALSE) +
  geom_point(aes(color = class), size = 4, show.legend = FALSE) +
  scale_color_manual(values = class_color) +
  base_theme +
  geom_text(aes(x = number, y = class, label = number), nudge_y = -0.2) +
  labs(x = "Feature number", y = "")
plot

ggsave(plot,
       filename = "omics_feature_number.pdf",
       width = 7,
       height = 7)
