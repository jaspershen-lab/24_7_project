no_function()

library(tidyverse)
setwd(masstools::get_project_wd())
rm(list = ls())

source("1_code/tools.R")
source("1_code/laggedCor/1-calculate_lagged_correlation.R")
source("1_code/laggedCor/2-lagged_cor_result.R")
source("1_code/laggedCor/3-evaluate_lagged_cor.R")
source("1_code/laggedCor/4-extract_data.R")
source("1_code/laggedCor/5-lagged_alignment_plot.R")
source("1_code/laggedCor/6-lagged_scatter_plot.R")
source("1_code/laggedCor/7-sample_matching_plot.R")
source("1_code/laggedCor/8-time_plot.R")
source("1_code/laggedCor/9-utils.R")

###data loading
####load all omics data loess data

####load all omics data loess data
load(
  here::here(
    "3_data_analysis/combine_omics/data_preparation/new_expression_data"
  )
)
load(here::here(
  "3_data_analysis/combine_omics/data_preparation/new_sample_info"
))
load(here::here(
  "3_data_analysis/combine_omics/data_preparation/new_variable_info"
))

variable_info <-
  new_variable_info

expression_data <-
  new_expression_data

sample_info =
  new_sample_info

rownames(expression_data) == variable_info$variable_id
colnames(expression_data) == sample_info$sample_id

######
dir.create("3_data_analysis/inter_omics_correlation")
setwd("3_data_analysis/inter_omics_correlation/")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between metabolic_panels
#global correlation

###get the matched index
time1 = sample_info$accurate_time
time2 = sample_info$accurate_time
x = as.numeric(expression_data[1,])
y = as.numeric(expression_data[1,])
time_tol = 300 / 60 # the biggest time shift is 5 hours
step = 30 / 60 # step is 0.5 hour

all_idx <-
  calculate_match_index(
    x = x,
    y = y,
    time1 = time1,
    time2 = time2,
    time_tol = time_tol,
    step = step,
    threads = 10
  )

all_idx

plot <-
  match_index_plot1(object = all_idx)
plot

# ggsave(plot, filename = "match_plot1.pdf", width = 9, height = 7)
# ggsave(plot, filename = "match_plot1.png", width = 9, height = 7)

plot =
  match_index_plot2(object = all_idx)
plot
# ggsave(plot, filename = "match_plot2.pdf", width = 9, height = 7)
# ggsave(plot, filename = "match_plot2.png", width = 9, height = 7)

dir.create("lagged_correlation")
total_number = nrow(expression_data)

library(future)
plan(multisession, workers = 2)
library(magrittr)

# index_list <- 
#   list(c(1:100),
#        c(101:300),
#        c(301:600),
#        c(601:1000),
#        c(1001:(total_number - 1)))  
# 
# 1:length(index_list) %>%
#   purrr::map(function(index) {
#     lagged_result <-
#       purrr::map(index_list[[index]], function(i) {
#         cat(i, " ")
#         temp_result <-
#           purrr::map((i + 1):total_number,
#                      .f = function(j) {
#                        cat(paste(i, j, sep = "-"), " ")
#                        x = as.numeric(expression_data[i,])
#                        time1 = sample_info$accurate_time
#                        y = as.numeric(expression_data[j,])
#                        time2 = sample_info$accurate_time
#                        
#                        result <-
#                          calculate_lagged_correlation(
#                            x = x,
#                            y = y,
#                            time1 = time1,
#                            time2 = time2,
#                            time_tol = 300 / 60,
#                            step = 30 / 60,
#                            min_matched_sample = 10,
#                            all_idx = all_idx$idx,
#                            progressbar = FALSE
#                          )
#                        
#                        p = result@all_cor_p
#                        p = p.adjust(p,
#                                     method = "BH",
#                                     n = nrow(variable_info) * (nrow(variable_info) - 1) / 2)
#                        if (all(p > 0.05, na.rm = TRUE)) {
#                          result = NULL
#                        }
#                        result
#                      }
#           )
#         names(temp_result) = rownames(expression_data)[(i + 1):total_number]
#         temp_result
#       })
#     
#     names(lagged_result) <-
#       rownames(expression_data)[index_list[[index]]]
#     save(lagged_result, file = file.path("lagged_correlation",
#                                          paste("lagged_result", index, sep = "_")))
#   })



####load lagged cor results
load("lagged_correlation/lagged_result_1")

# #####evaluate peak quality
#

evaluate_lagged_cor(object = result)

evaluate_lagged_cor(object = lagged_result[[2]][[15]])

###output the matched plot
dir.create("sample_matching_plot")

idx =
  lagged_result[[2]] %>% lapply(class) %>% unlist() %>%
  `==`("list") %>%
  which()


load(here::here("3_data_analysis/summary_info/day_night_df"))


# cor_data =
#   purrr::map(1:length(lagged_result), function(i) {
#     cat(i, " ")
#     x = lagged_result[[i]]
#     purrr::map(1:length(x), function(j) {
#       # cat(j, " ")
#       y = x[[j]]
#       if (is.null(y) | length(y) == 1) {
#         return(NULL)
#       }
#       temp_data =
#       data.frame(
#         from = names(lagged_result)[i],
#         to = names(x)[j],
#         shift_time = y$shift_time,
#         cor = y$all_cor,
#         p = y$all_cor_p
#       )
#
#       result = evaluate_peak_quality(object = y, plot = FALSE)
#       temp_data$score = result$score
#       temp_data
#     }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(cor_data, file = "lagged_correlation/cor_data", compress = "xz")
load("lagged_correlation/cor_data")

##output the shift time vs correlation plots
temp_data =
  cor_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::group_by(from, to) %>%
  dplyr::filter(abs(cor) == max(abs(cor))) %>%
  dplyr::ungroup()

temp_data$p_adjust = p.adjust(temp_data$p,
                              method = "bonferroni",
                              n = nrow(expression_data) * (nrow(expression_data) - 1) /
                                2)
temp_data =
  temp_data %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::filter(shift_time != "(-15,15]")

dir.create("shift_time_vs_cor")

# for(i in 1:nrow(temp_data)){
#   cat(i, " ")
#   from = temp_data$from[i]
#   to = temp_data$to[i]
#   idx1 = which(names(lagged_result) == from)
#   idx2 = which(names(lagged_result[[idx1]]) == to)
#   result =
#     evaluate_peak_quality(object = lagged_result[[idx1]][[idx2]], plot = TRUE)
#   name =
#     paste(round(result$score, 4), "_",from, "_",to, ".pdf",sep = "") %>%
#     stringr::str_replace_all("/", "_")
#   ggsave(result$plot, filename = file.path("shift_time_vs_cor", name), width = 9, height = 7)
# }

###check the shift time vs correlation, and set the score cutoff as 0.5

lagged_cor =
  cor_data %>%
  dplyr::filter(!is.na(cor))

lagged_cor$cor[lagged_cor$shift_time != "(-15,15]" &
                 lagged_cor$score < 0.5] = 0

lagged_cor =
  lagged_cor %>%
  dplyr::group_by(from, to) %>%
  dplyr::filter(abs(cor) == max(abs(cor))) %>%
  dplyr::ungroup()

lagged_cor$p_adjust = p.adjust(lagged_cor$p,
                               method = "bonferroni",
                               n = nrow(expression_data) * (nrow(expression_data) - 1) /
                                 2)

lagged_cor =
  lagged_cor %>%
  dplyr::filter(p_adjust < 0.05)

global_cor =
  cor_data %>%
  dplyr::filter(shift_time == "(-15,15]")

global_cor$p_adjust = p.adjust(global_cor$p, method = "bonferroni",
                               nrow(expression_data) * (nrow(expression_data) - 1) /
                                 2)

global_cor =
  global_cor %>%
  dplyr::filter(p_adjust < 0.05)

dim(lagged_cor)
dim(global_cor)

table(lagged_cor$shift_time)

lagged_cor =
  lagged_cor %>%
  dplyr::mutate(name = paste(from, to, sep = "_")) %>%
  dplyr::select(name, dplyr::everything()) %>%
  dplyr::arrange(name)

global_cor =
  global_cor %>%
  dplyr::mutate(name = paste(from, to, sep = "_")) %>%
  dplyr::select(name, dplyr::everything()) %>%
  dplyr::arrange(name)

global_cor %>%
  dplyr::arrange(abs(cor)) %>%
  tail()

plot(as.numeric(expression_data["HILIC_POS_10.33_172.1332m/z",]),
     as.numeric(expression_data["HILIC_POS_11.00_87.0441m/z",]))

dim(global_cor)
dim(lagged_cor)

library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "lagged correlation",
#              gridLines = TRUE)
# addWorksheet(wb, sheetName = "global correlation",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = lagged_cor,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 2, x = global_cor,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "lagged_correlation/cor_data.xlsx", overwrite = TRUE)
