no_function()
library(tidyverse)
setwd(masstools::get_project_wd())
rm(list = ls())
source("1_code/tools.R")
# source("1_code/modified_dtw.R")
source("1_code/lagged_correlation.R")

{
  load(here::here("3_data_analysis/summary_info/day_night_df"))
  
  ####load data (combined omics data)
  load("3_data_analysis/combine_omics/data_preparation/new_sample_info")
  load("3_data_analysis/combine_omics/data_preparation/new_variable_info")
  load("3_data_analysis/combine_omics/data_preparation/new_expression_data")
  
  expression_data = new_expression_data
  sample_info = new_sample_info
  variable_info = new_variable_info
  
  load("3_data_analysis/combine_omics/data_preparation/cortisol_variable_info")
  load("3_data_analysis/combine_omics/data_preparation/cytokine_variable_info")
  load("3_data_analysis/combine_omics/data_preparation/lipidomics_variable_info")
  load("3_data_analysis/combine_omics/data_preparation/metabolomics_variable_info")
  load("3_data_analysis/combine_omics/data_preparation/proteomics_variable_info")
  load("3_data_analysis/combine_omics/data_preparation/total_protein_variable_info")
  load(
    "3_data_analysis/combine_omics/data_preparation/metabolic_panel_variable_info"
  )
}

dir.create("3_data_analysis/k_means_clustering", recursive = TRUE)
setwd("3_data_analysis/k_means_clustering")

table(variable_info$data_type)

dim(expression_data)

sample_info %>%
  ggplot(aes(accurate_time, time)) +
  geom_point() +
  base_theme

library(Mfuzz)

temp_data <-
  expression_data

rownames(temp_data)

time <- c(1:ncol(temp_data))

temp_data <- rbind(time, temp_data)

temp_data2 =
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

row.names(temp_data)[1] <- "time"

# write.table(
#   temp_data,
#   file = "temp_data.txt",
#   sep = '\t',
#   quote = FALSE,
#   col.names = NA
# )

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

# plot <-
# Dmin(
#   data.s,
#   m = m1,
#   crange = seq(2, 40, 1),
#   repeats = 3,
#   visu = TRUE
# )
#
# plot <-
# plot %>%
#   data.frame(distance = plot,
#              k = seq(2,40,1)) %>%
#   ggplot(aes(k, distance)) +
#   geom_point(shape = 21, size = 4, fill = "black") +
#   geom_smooth() +
#   geom_segment(aes(x = k, y = 0, xend = k, yend = distance)) +
#   theme_bw() +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) +
#   labs(
#     x = "Cluster number",
#     y = "Min. centroid distance"
#   ) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#
# plot
#
# ggsave(plot, filename = "distance_k_number.pdf", width = 7, height = 7)

clust = 11

# c <- mfuzz(data.s, c = clust, m = m1)
#
# mfuzz.plot(eset = data.s,
#            min.mem = 0.8,
#            cl = c,
#            mfrow=c(3,4),
#            time.labels = time,
#            new.window = FALSE)
#
# names(c$cluster) <- rownames(temp_data2)[-1]
# rownames(c$membership) <- rownames(temp_data2)[-1]
# save(c, file = "c")
load("c")

####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
center <- c$centers
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

acore <- acore(data.s, c, min.acore = 0)
acore

centers <- c$centers
names(c$cluster) == rownames(c$membership)

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

# openxlsx::write.xlsx(x = cluster_info,
#                      file = "cluster_info.xlsx", asTable = TRUE)

#####output the expression data of different clusters

###plot for each cluster
# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
#   dir.create(paste("cluster", cluster_idx, sep = "_"))
#   cluster_data <-
#     cluster_info %>%
#     dplyr::filter(cluster_idx == cluster_idx) %>%
#     dplyr::select(1, 1 + cluster_idx)
#
#   colnames(cluster_data) <- c("variable_id", "membership")
#
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.9)
#
#   openxlsx::write.xlsx(
#     x = cluster_data,
#     file = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".xlsx", sep = "")
#     ),
#     asTable = TRUE,
#     overwrite = TRUE
#   )
#
#   ###cluster plot
#
#   temp =
#     temp_data2[cluster_data$variable_id, ] %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership),
#       names_to = "sample_id",
#       values_to = "value"
#     ) %>%
#     dplyr::left_join(sample_info[, c("sample_id", "accurate_time")], by = "sample_id")
#
#   plot <-
#     ggplot() +
#     geom_rect(
#       mapping = aes(
#         xmin = start,
#         xmax = end,
#         ymin = -Inf,
#         ymax = Inf
#       ),
#       fill = "lightyellow",
#       data = day_night_df %>%
#         dplyr::filter(as.character(day) %in% sample_info$day),
#       show.legend = FALSE
#     ) +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(accurate_time, value, group = variable_id, color = membership),
#               data = temp) +
#     scale_x_datetime(
#       breaks = scales::date_breaks("12 hour"),
#       date_labels = "%a %H:%M",
#       timezone = "America/Los_Angeles"
#     ) +
#     scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1,
#         vjust = 1,
#         size = 12
#       ),
#       panel.background = element_rect(fill = alpha("grey", 0.2)),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste(
#         "Cluster ",
#         cluster_idx,
#         " (",
#         nrow(cluster_data),
#         " molecules)",
#         sep = ""
#       )
#     )
#
#   plot
#
#   ggsave(
#     plot,
#     filename = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".pdf", sep = "")
#     ),
#     width = 20,
#     height = 7
#   )
# }

## (4) feature number
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx")
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx")
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx")
cluster4 <- readxl::read_xlsx("cluster_4/cluster4.xlsx")
cluster5 <- readxl::read_xlsx("cluster_5/cluster5.xlsx")
cluster6 <- readxl::read_xlsx("cluster_6/cluster6.xlsx")
cluster7 <- readxl::read_xlsx("cluster_7/cluster7.xlsx")
cluster8 <- readxl::read_xlsx("cluster_8/cluster8.xlsx")
cluster9 <- readxl::read_xlsx("cluster_9/cluster9.xlsx")
cluster10 <- readxl::read_xlsx("cluster_10/cluster10.xlsx")
cluster11 <- readxl::read_xlsx("cluster_11/cluster11.xlsx")

###annotation for each cluster
cluster1 = data.frame(cluster1, cluster = "1")
cluster2 = data.frame(cluster2, cluster = "2")
cluster3 = data.frame(cluster3, cluster = "3")
cluster4 = data.frame(cluster4, cluster = "4")
cluster5 = data.frame(cluster5, cluster = "5")
cluster6 = data.frame(cluster6, cluster = "6")
cluster7 = data.frame(cluster7, cluster = "7")
cluster8 = data.frame(cluster8, cluster = "8")
cluster9 = data.frame(cluster9, cluster = "9")
cluster10 = data.frame(cluster10, cluster = "10")
cluster11 = data.frame(cluster11, cluster = "11")

cluster =
  rbind(
    cluster1,
    cluster2,
    cluster3,
    cluster4,
    cluster5,
    cluster6,
    cluster7,
    cluster8,
    cluster9,
    cluster10,
    cluster11
  )

variable_info =
  variable_info %>%
  dplyr::left_join(cluster, by = "variable_id")

variable_info$cluster[is.na(variable_info$cluster)] = "Other"

variable_info$mol_name[!is.na(variable_info$Lipid_Name)] =
  variable_info$Lipid_Name[!is.na(variable_info$Lipid_Name)]

table(variable_info$cluster)

####Heatmap to show the clusters
library(ComplexHeatmap)
dim(expression_data)
rownames(expression_data) == cluster_info$variable_id

temp_data = expression_data[cluster_info$variable_id, ]

temp_cluster_info =
  cluster_info[match(rownames(temp_data), cluster_info$variable_id),]

rownames(expression_data) == temp_cluster_info$variable_id

cluster = temp_cluster_info$cluster

library(lubridate)

temp_data =
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>%
  t() %>%
  as.data.frame()

range(temp_data, na.rm = TRUE)
temp_data[temp_data > 3] = 3
temp_data[temp_data < -3] = -3

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

table(cluster)

library(circlize)

col_fun = colorRamp2(
  breaks = c(-3, 0, 3),
  colors =
    c("#366A9FFF", "white", "red"),
  transparency = 0
)

table(cluster)

plot =
  Heatmap(
    temp_data,
    col = col_fun,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    name = "Z-score",
    border = TRUE,
    column_names_gp = gpar(fontsize = 10, angle = 45),
    column_names_rot = 45,
    column_labels = x_labels,
    row_split = factor(cluster, levels = unique(cluster)),
    # row_title = rep("", length(unique(cluster))),
    row_title_rot = 0,
    na_col = "grey"
    # right_annotation = rowAnnotation(foo = anno_block(
    #   labels = unique(cluster),
    #   labels_gp = gpar(
    #     col = "black", fontsize = 10
    #   )))
  )
rownames(temp_data)[unlist(row_order(plot))]
rownames(temp_data)[unlist(row_order(plot))] == temp_cluster_info$variable_id

# plot

temp_idx =
  row_order(plot) %>%
  purrr::map(function(x) {
    rownames(temp_data)[x]
  })

temp_data =
  temp_data[unlist(temp_idx),]

#####use ggplot2 to get the heatmap
colnames(temp_data)

# plot1 = ggplotify::as.ggplot(plot)
#
# plot1
#
# ggsave(plot1, filename = "heatmap.dent.pdf", width = 7, height = 7)

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
  dplyr::left_join(cluster_info[, c("variable_id", "cluster")],
                   by = c("variable_name" = "variable_id")) %>%
  tibble::rowid_to_column(var = "variable_id") %>%
  dplyr::select(-variable_name) %>%
  dplyr::select(variable_id, cluster, everything()) %>%
  dplyr::mutate(variable_id = as.numeric(variable_id)) %>%
  tidyr::pivot_longer(
    cols = -c(variable_id, cluster),
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
      cluster = temp_data1_na$cluster[idx],
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

plot <-
  ggplot(temp_data1,
         aes(
           xmin = time1,
           xmax = time2,
           ymin = variable_id,
           ymax = variable_id + 1
         )) +
  geom_rect(aes(fill = fill), colour = NA) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
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
  facet_grid(
    rows = vars(cluster),
    cols = vars(date),
    scales = "free",
    space = "free"
  )
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
