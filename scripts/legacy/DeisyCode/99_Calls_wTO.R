args=(commandArgs(TRUE))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
print(N)
print(save)
print(file)

source("./99_newwTO.R")
input = fread(file) %>%
  as.data.frame()

wto = wTO.faster(Data = input, n = N)

wto %>% fwrite(save)

message(paste("Saved calculations to\n", save))
