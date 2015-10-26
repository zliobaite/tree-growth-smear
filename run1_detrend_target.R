# 2015 03 19 I.Zliobaite
# zliobaite@gmail.com
# constructs response variables out of individual trees

library(dplR)

input_file_raw_widths <- 'input_data/raw_widths.csv'
input_file_raw_heights <- 'input_data/raw_heights.csv' #height growth

output_file_target_widths <- 'processed_data/targets_widths.csv'
output_file_target_heights <- 'processed_data/targets_heights.csv'

detrend_method_widths <- "ModNegExp"
detrend_method_heights <- "Spline"

data_widths <- read.csv(input_file_raw_widths,sep = ',',header = FALSE)
data_heights <- read.csv(input_file_raw_heights,sep = ',',header = FALSE)

year_start_detrend_widths = 1984
year_end_detrend_widths = 2014
year_start_target_widths = 1998
year_end_target_widths = 2013

year_start_detrend_heights = 1982
year_end_detrend_heights = 2012
year_start_target_heights = 1998
year_end_target_heights = 2012

c = 6; #for biwariate mean
cik = 10 #for biwariate mean
#--------------------------------------------------------
#read data
data_widths <- read.csv(input_file_raw_widths,sep = ',',header = FALSE)
data_heights <- read.csv(input_file_raw_heights,sep = ',',header = FALSE)
#--------------------------------------------------------
# detrend data
detrend_data <- function(data0,year_start,year_end,detrend_method)
{
  p <- dim(data0)[2]
  data0 <- data.frame(data0)
  row.names(data0) <- data0[,1]
  data0 <- data0[as.character(c(year_start:year_end)),c(2:p)] #select relevant
  class(data0) <- c("rwl", "data.frame")
  data_detrended <- detrend(rwl = data0, method = detrend_method, make.plot = TRUE)
  return(data_detrended)
}
#--------------------------------------------------------
widths_detrended <- detrend_data(data_widths,year_start_detrend_widths,year_end_detrend_widths,detrend_method_widths)

heights_detrended <- detrend_data(data_heights,year_start_detrend_heights,year_end_detrend_heights,detrend_method_heights)

#--------------------------------------------------------
estimate_bi_mean <- function(data,c,cik)
{
  ystar <- apply(data,1,mean)
  n_trees <- dim(data)[2]
  estimates <- ystar
  for (sk in 1:cik)
  {
    s <- t(data) - matrix(1,n_trees,1)%*%ystar
    s <- apply(abs(s),2,median)
    b <- (t(data) - matrix(1,n_trees,1)%*%ystar)/(c*matrix(1,n_trees,1)%*%s)
    b <- b^2
    w <- (1 - b)^2
    ind <- which(w>1)
    w[ind] <- 0
    ystar <- apply(w*t(data),2,sum)/apply(w,2,sum)
  }
  return(ystar)
}

make_target <- function(data_detrended,year_start_target,year_end_target,c,cik,output_file_target)
{
  years_target <- c(year_start_target:year_end_target)
  data_detrended <- data_detrended[as.character(years_target),]
  bi_mean <- estimate_bi_mean(data_detrended,c,cik)
  #detrend linear
  fit <- lm(bi_mean ~ years_target)
  trend <- coef(fit)[2]*years_target + coef(fit)[1]
  bi_mean <- bi_mean - trend
  bi_mean <- bi_mean / sd(bi_mean)
  target <- cbind(years_target,bi_mean)
  target <- round(target,digits = 4)
  write.table(target, file = output_file_target, row.names=FALSE, col.names=FALSE, sep="\t")
}

make_target(widths_detrended,year_start_target_widths,year_end_target_widths,c,cik,output_file_target_widths)

make_target(heights_detrended,year_start_target_heights,year_end_target_heights,c,cik,output_file_target_heights)
