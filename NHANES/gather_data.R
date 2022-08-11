library(tidyverse)
xpt_files  <- list.files(pattern="*\\.xpt", ignore.case=TRUE)

df_list  <- lapply(xpt_files, function(x) haven::read_xpt(x))
names(df_list) <- xpt_files

df <- df_list %>% reduce(left_join, by="SEQN")

write_csv(df, file="nhanes_data.csv")