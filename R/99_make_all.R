#' ---
#' title: Make All
#' author: Keith Baggerly
#' date: "`r Sys.Date()`"
#' output: github_document
#' ---

library(here)
library(rmarkdown)
library(fs)
library(tools)

if(!dir.exists(here::here("results"))){
  dir.create(here::here("results"))
}

files_in_r_to_run <- 
  c("01_gather_raw_data.R",
    "02_parse_potti_data.Rmd",
    "03_parse_nci60_data.Rmd",
    "04_parse_geo_data.Rmd",
    "05_identify_potti_sources.Rmd",
    "06_report_potti_data_sources.Rmd",
    "90_kludges_and_warnings.Rmd")


#################
## Generate new files
#################

file_prefixes <- 
  tools::file_path_sans_ext(files_in_r_to_run)

for(i1 in 1:length(files_in_r_to_run)){
    
  rmarkdown::render(
    here::here("R", files_in_r_to_run[i1]),
    output_format = 
      github_document(html_preview = TRUE, toc = TRUE),
    envir = new.env())
  
  files_to_move <- basename(dir_info(here::here("R"))$path)
  files_to_move <- 
    files_to_move[grep(file_prefixes[i1], files_to_move)]
  files_to_move <- 
    files_to_move[-which(files_to_move == files_in_r_to_run[i1])]
  cached_file_index <- grep("cache$", files_to_move)
  if(length(cached_file_index) > 0){
    files_to_move <- 
      files_to_move[-cached_file_index]
  }
  
  files_to_overwrite <- 
    here::here("results", files_to_move)
  files_to_overwrite <- 
    files_to_overwrite[file_exists(files_to_overwrite)]
  
  if(length(files_to_overwrite) > 0){
    fs::file_delete(files_to_overwrite)
  }

  fs::file_move(here::here("R", files_to_move), 
                here::here("results"))
  
}

rmarkdown::render(
  here::here("README.Rmd"),
  output_format = 
    github_document(html_preview = TRUE, toc = TRUE),
  envir = new.env())

## Do NOT use the output_dir option to rmarkdown::render,
## since this replaces relative paths with absolute ones.
