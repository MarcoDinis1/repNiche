library("ecospat")

lst <- list.files(".", pattern = "\\.csv$", full.names = T) # name long/lat columns as x/y beforehand and put them in positions 3:4
lst2 <- list.files(".", pattern = "\\.csv$")


for (i in 1:length(lst))
{
  input <- read.csv(lst[i], sep = ";")
  input$x <- as.numeric(input$x)
  input$y <- as.numeric(input$y)
  input2 <- unique(input[,c(3:4)])
  output <- ecospat.occ.desaggregation(input2, min.dist = 0.041666667) # set min.dist to 0.0083333333 for 1km pixel and 0.041666667 for 4km
  write.csv(output, file = paste(".\\", gsub(".csv", "_dec4km", lst2[i]), ".csv", sep = ""))
}


