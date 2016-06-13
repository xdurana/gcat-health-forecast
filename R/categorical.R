library(xlsx)

dictionary <- read.xlsx2('output/GCAT Health Forecast.data-dictionary.xlsx', sheetIndex = 2)
categories <- as.character(unique(dictionary$variable))
write.table(categories, 'config/categorical.txt', row.names = FALSE, quote = FALSE)
