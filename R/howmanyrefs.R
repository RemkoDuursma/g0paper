

rl <- readLines("manuscript.Rmd")

library(stringr)
ref <- str_extract_all(rl, "@[a-zA-Z]+[0-9]{4}")
length(unique(do.call(c, ref)))
