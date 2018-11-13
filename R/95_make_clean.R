#' ---
#' title: Make Clean
#' author: Keith Baggerly
#' date: "`r Sys.Date()`"
#' output: github_document
#' ---

library(here)
library(fs)

dirs_to_clean <- c("results")

for(i1 in 1:length(dirs_to_clean)){
  temp_file_list <- 
    dir(here::here(dirs_to_clean[i1]))
  file_delete(here::here(dirs_to_clean[i1],temp_file_list))
}
