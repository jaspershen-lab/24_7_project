no_function()
library(tidyverse)
setwd(masstools::get_project_wd())
rm(list = ls())
source("1_code/tools.R")

{
  load(here::here("3_data_analysis/summary_info/day_night_df"))
  
  ####load data (combined omics data)
  load("3_data_analysis/food/data_preparation/food_sample_info")
  load("3_data_analysis/food/data_preparation/food_expression_data")
  load("3_data_analysis/food/data_preparation/food_variable_info")
  
  expression_data = food_expression_data
  sample_info = food_sample_info
  variable_info = food_variable_info
}

dir.create("3_data_analysis/food/data_overview")
setwd("3_data_analysis/food/data_overview")

####change some units
variable_info$mol_name

expression_data[8,] <-
  expression_data[8,] / 1000

expression_data[9,] <-
  expression_data[9,] / 1000

expression_data[12,] <-
  expression_data[12,] / 1000

variable_info$mol_name[c(8, 9, 12)] <-
  variable_info$mol_name[c(8, 9, 12)] %>%
  stringr::str_replace("mg", "g")

####Heatmap to show the clusters
library(ComplexHeatmap)
dim(expression_data)
rownames(expression_data) == variable_info$variable_id

temp_data = expression_data

rownames(temp_data) <-
  variable_info$mol_name

library(lubridate)

temp_data =
  temp_data

###for each nutriton, scale them to 0 - 100
temp_data2 <-
  temp_data %>%
  apply(1, function(x) {
    x <- as.numeric(x)
    (x - min(x)) / (max(x) - min(x)) * 1
  }) %>%
  t() %>%
  as.data.frame()

colnames(temp_data2) <-
  colnames(temp_data)

temp_data <-
  temp_data2

range(temp_data, na.rm = TRUE)

dim(temp_data)

dim(sample_info)

x_labels =
  colnames(temp_data) %>%
  stringr::str_replace("\\:00", "") %>%
  stringr::str_replace("2019-04-29", "Mon") %>%
  stringr::str_replace("2019-04-30", "Tue") %>%
  stringr::str_replace("2019-05-01", "Wed") %>%
  stringr::str_replace("2019-05-02", "Thu") %>%
  stringr::str_replace("2019-05-03", "Fri") %>%
  stringr::str_replace("2019-05-04", "Sat") %>%
  stringr::str_replace("2019-05-05", "Sun") %>%
  stringr::str_replace("2019-05-06", "Mon") %>%
  stringr::str_replace("2019-05-07", "Tue")

x_labels[-seq(1, length(x_labels), by = 5)] = ""

library(circlize)

range(temp_data)

col_fun = colorRamp2(
  breaks = seq(range(temp_data)[1], range(temp_data)[2], length.out = 9),
  colors =
    RColorBrewer::brewer.pal(n = 9, name = "Blues"),
  transparency = 0
)

plot =
  Heatmap(
    temp_data,
    col = col_fun,
    show_row_names = TRUE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    name = "Food",
    border = TRUE,
    column_names_gp = gpar(fontsize = 10, angle = 45),
    column_names_rot = 45,
    column_labels = x_labels,
    # row_title = rep("", length(unique(cluster))),
    row_title_rot = 0,
    na_col = "grey",
    rect_gp = gpar(col = "grey"),
    # right_annotation = rowAnnotation(foo = anno_block(
    #   labels = unique(cluster),
    #   labels_gp = gpar(
    #     col = "black", fontsize = 10
    #   )))
  )

plot

#####use ggplot2 to get the heatmap
colnames(temp_data)

plot1 = ggplotify::as.ggplot(plot)

plot1

ggsave(plot1,
       filename = "heatmap.dent.pdf",
       width = 7,
       height = 1.5)

ggsave(plot1,
       filename = "heatmap.dent.png",
       width = 7,
       height = 1.5)

# time = lubridate::as_datetime(colnames(temp_data), tz = "America/Los_Angeles")
time = sample_info$accurate_time

min_time = lubridate::as_datetime("2019-04-29 00:00:01", tz = "America/Los_Angeles")
max_time = lubridate::as_datetime("2019-05-07 23:59:00", tz = "America/Los_Angeles")

# time[length(time)] = time[length(time)] + 30*60

colnames(temp_data)

time1 = c(time[1] - 30 * 60, time[-length(time)])
time2 = c(time)

names(time1) =
  names(time2) =
  colnames(temp_data)

time =
  data.frame(sample_name = colnames(temp_data),
             time1,
             time2)

####for the time when there are no sampling, just set it as NA
time =
  time %>%
  dplyr::mutate(width = as.numeric(difftime(time2, time1, units = "min")))

library(plyr)

time =
  time %>%
  plyr::dlply(.variables = .(sample_name)) %>%
  purrr::map(function(x) {
    if (x$width == 30) {
      return(x)
    } else{
      x1 = x
      x2 = x
      x1$time2 = x1$time1 + 30 * 60
      x1$width = as.numeric(difftime(x1$time2, x1$time1, units = "min"))
      x2$time1 = x1$time2
      x2$width = as.numeric(difftime(x2$time2, x2$time1, units = "min"))
      x2$sample_name = paste(x2$sample_name, 1, sep = "_")
      rbind(x1, x2)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data1 =
  temp_data

add_data = matrix(NA,
                  nrow = nrow(temp_data1),
                  ncol = length(time$sample_name[which(time$width > 30)])) %>%
  as.data.frame()

colnames(add_data) = time$sample_name[which(time$width > 30)]

temp_data1 =
  cbind(temp_data1, add_data)

temp_data1 =
  temp_data1 %>%
  tibble::rownames_to_column(var = "variable_name") %>%
  tibble::rowid_to_column(var = "variable_id") %>%
  dplyr::select(-variable_name) %>%
  dplyr::select(variable_id, everything()) %>%
  dplyr::mutate(variable_id = as.numeric(variable_id)) %>%
  tidyr::pivot_longer(
    cols = -c(variable_id),
    names_to = "sample_name",
    values_to = "fill"
  ) %>%
  dplyr::left_join(time, by = "sample_name")

sum(is.na(temp_data1$fill))

range(temp_data1$fill, na.rm = TRUE)

####add more data
temp_data1_normal <-
  temp_data1 %>%
  dplyr::filter(width == 30)

temp_data1_na <-
  temp_data1 %>%
  dplyr::filter(width > 30)

temp_data1_na <-
  seq_len(nrow(temp_data1_na)) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    new_time <-
      seq(temp_data1_na$time1[idx], temp_data1_na$time2[idx], by = 30 *
            60)
    data.frame(
      variable_id = temp_data1_na$variable_id[idx],
      sample_name = temp_data1_na$sample_name[idx],
      fill = temp_data1_na$fill[idx],
      time1 = new_time[-length(new_time)],
      time2 = new_time[-1],
      width = 30
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data1 <-
  rbind(temp_data1_normal,
        temp_data1_na)

date <-
  format(temp_data1$time2, "%Y-%m-%d")

temp_data1$date <- date

temp_data1[is.na(temp_data1)] <- 0

plot <-
  ggplot(temp_data1,
         aes(
           xmin = time1,
           xmax = time2,
           ymin = variable_id,
           ymax = variable_id + 1
         )) +
  geom_rect(aes(fill = fill), colour = NA) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
                       na.value = "grey") +
  scale_x_datetime(
    breaks = scales::date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    expand = expansion(mult = c(0, 0)),
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 10
    ),
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_rect(fill = alpha("grey", 0.2)),
    panel.spacing = unit(0.1, "lines")
  ) +
  facet_grid(cols = vars(date),
             scales = "free",
             space = "free")
plot
# ggsave(plot,
#        filename = "molecular_heatmap2.pdf",
#        width = 10,
#        height = 7)

#####
####mosaic plot to show the difference
library(ggmosaic)

temp_data2 <-
  variable_info %>%
  dplyr::filter(cluster != "Other") %>%
  dplyr::mutate(cluster = factor(cluster, levels = c(as.character(c(
    1:11
  )))))

plot2 =
  temp_data2 %>%
  ggplot() +
  geom_mosaic(aes(x = ggmosaic::product(data_type, cluster), fill = data_type),
              offset = 0.01) +
  theme_mosaic() +
  scale_fill_manual(values = class_color) +
  labs(x = "", y = "") +
  theme(panel.border = element_rect(color = "black",
                                    fill = "transparent"))

plot2

# ggsave(plot2, filename = "cluster_data_type.pdf", width = 7, height = 3)
