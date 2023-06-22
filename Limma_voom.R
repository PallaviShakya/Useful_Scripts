## Set working directory
setwd("~/git_repositories/tag_seq/star_align_all")

rm(list = ls())
options(scipen = 999)

install.packages("devtools")

install.packages("glue")
install.packages("tidyverse")

library("tidyverse")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "tidyverse", "stephenturner/annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx", "rstudio/gt", "plyr", "glue", "Glimma", "sva"))
stopifnot(suppressMessages(sapply(c("edgeR", "tidyverse", "annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx", "gt", "plyr", "glue", "Glimma", "sva"),
                                  require, character.only = TRUE)))

sampleNames <- list.files(path = glue :: glue(getwd()), pattern = "*ReadsPerGene.out.tab")%>%
  stringr::str_split_fixed("_", n = 4) %>%
  tibble::as_tibble() %>%
  tidyr::unite(Name, c(V1:V3), sep = "_") %>%
  dplyr::select(Name) %>% 
  purrr::flatten_chr()

