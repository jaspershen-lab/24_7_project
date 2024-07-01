##
no_function()
setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())
source(here::here("1_code/tools.R"))

dir.create("3_data_analysis/food/data_preparation", recursive = TRUE)
setwd("3_data_analysis/food/data_preparation")

load("expression_data")
load("sample_info")
load("variable_info")

dim(expression_data)
dim(sample_info)
dim(variable_info)

head(sample_info)
head(variable_info)

sum(is.na(expression_data))

which(is.na(expression_data), arr.ind = TRUE)

expression_data[which(is.na(expression_data), arr.ind = TRUE)] <-0

food_expression_data <- 
  expression_data

food_sample_info <- 
  sample_info %>% 
  dplyr::mutate(class = "Subject")

food_variable_info <- 
  variable_info

library(massdataset)

food_data <-
create_mass_dataset(expression_data = food_expression_data,
                    sample_info = food_sample_info,
                    variable_info = food_variable_info)

base::save(food_data, file = "food_data")
base::save(food_expression_data, file = "food_expression_data")
base::save(food_sample_info, file = "food_sample_info")
base::save(food_variable_info, file = "food_variable_info")



