

url <- "http://onlinelibrary.wiley.com/doi/10.1111/pce.12871/full"

library(rvest)
library(stringr)

tab <- read_html(url, encoding="latin1")

tab <- html_nodes(tab, "table")

table1 <- html_table(tab[1])

df <- as.data.frame(table1[[1]][-1,]) 
names(df) <- c("species","growthcond","model","m","g0","method","R2","ref")

df$m <- gsub("\\(.+\\)","", df$m)

df$m <- sapply(strsplit(df$m, "±"), "[", 1)

df$m <- gsub(" a", "", df$m)
df$m <- gsub("a", "", df$m)
df$m <- gsub("\U00002009", "-", df$m)
df$m <- gsub("≈", "", df$m)
df$m <- gsub(" to ", "-", df$m)

df$m <- str_trim(df$m)

df$m <- gsub("[^\\.0-9]", "--", df$m)
df$m <- gsub("^--", "", df$m)
df$m <- gsub("--$", "", df$m)

l <- strsplit(df$m, "--", fixed=TRUE)

df$m <- sapply(l, function(x){
  x <- as.numeric(x)
  if(length(x) == 1)
    return(x)
  else
    return((x[2] + x[1])/2)
})
df$m[which(df$m > 100)] <- df$m[which(df$m > 100)]/10

df <- subset(df, model == "BB")

g0 <- df$g0

g0 <- gsub("\\(.+\\)","", g0)

g0 <- gsub(" a", "", g0)
g0 <- gsub("a", "", g0)
g0 <- str_trim(g0)
g0 <- gsub("[^\\.0-9]", "--", g0)
l <- strsplit(g0, "\U00002009")
g0 <- sapply(l, "[", 1)


ii <- c(1,5,107,179,181) 
g0l <- strsplit(g0[ii], "[^\\.0-9]")

g0l <- sapply(g0l, function(x){
  x <- as.numeric(x)
  if(length(x) == 1)
    return(x)
  else
    return((x[3] + x[1])/2)
})
g0[ii] <- g0l

# the ones that don't convert; we don't want those!
df$g0 <- as.numeric(g0)
df$R2 <- as.numeric(df$R2)


with(df, plot(R2, g0))
l <- loess(g0 ~ R2, data=df)
plot_loess(l, add=TRUE)



