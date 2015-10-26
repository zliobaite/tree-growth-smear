# 2015 03 20 I.Zliobaite
# from 30 min SMEAR data makes 4 week windows

input_file_name <- 'input_data/raw_SMEAR.csv'
input_file_names_all <- 'input_data/names_SMEAR.csv'

output_file_features_all <- 'processed_data/data_aggregated_SMEAR.csv'
output_file_encoding_all <- 'processed_data/encoding_aggregated_SMEAR.csv'
output_file_year_range <- 'processed_data/year_aggregated_SMEAR.csv'

out_table1 <- 'results/table1.txt'
out_fet_select <- 'results/fet_select.csv'

param_aggregation_function <- mean

data_all <- read.csv(input_file_name, header = FALSE, sep = ',',na.string = 'NaN') #
#apply(data_all,2,max,na.rm=TRUE)



year_range <- unique(data_all[,1])
p <- dim(data_all)[2]

records_in_week <- 24*2*7
weeks_back <- 52
weeks_forward <- 36
step_weeks <- 2
window_weeks <- 4

vector_mean <- function(data1,param_aggregation)
{
  #print(dim(data1))
  vec_now <- apply(data1,2,param_aggregation_function,na.rm=TRUE)
  vec_na <- apply(is.na(data1),2,sum)
  ind_na <- which(vec_na>0)
  if (length(vec_na)>0)
  {
    vec_now[ind_na] <- NA
  }
  return(vec_now)
}

features_all <- c()
for (sk in 2:length(year_range))
{
  year_now <- year_range[sk]
  year_before <- year_range[sk-1]
  ind_year_now <- which(data_all[,1] == year_now) 
  ind_year_before <- which(data_all[,1] == year_before) 
  data_now <- data_all[ind_year_now,c(7:p)] #assume sorted
  data_before <- data_all[ind_year_before,c(7:p)] #assume sorted
  #init
  features_now <- c()
  encoding_fet <- c()
  encoding_week <- c()
  #year back
  time_no <- 1
  records_this_year <- dim(data_before)[1]
  ind_start <- records_this_year - weeks_back*records_in_week + 1
  ind_end <- ind_start + records_in_week*window_weeks - 1
  while (ind_end <= records_this_year)
  {    
    vec_now <- vector_mean(data_before[c(ind_start:ind_end),],param_aggregation_function)
    features_now <- c(features_now,vec_now)
    encoding_fet <- c(encoding_fet,c(1:(p-6)))
    encoding_week <- c(encoding_week,rep(1,(p-6))*time_no)
    time_no <- time_no + 1
    ind_start <- ind_start + records_in_week*step_weeks
    ind_end <- ind_end + records_in_week*step_weeks
  }
  # middle of year
  ind_end <- records_in_week*window_weeks/2
  ind_start <- records_this_year - records_in_week*window_weeks/2 + 1
  data_mid <- rbind(data_now[c(1:ind_end),],data_before[c(ind_start:records_this_year),])
  vec_now <- vector_mean(data_mid,param_aggregation_function)
  features_now <- c(features_now,vec_now)  
  encoding_fet <- c(encoding_fet,c(1:(p-6)))
  encoding_week <- c(encoding_week,rep(time_no,(p-6)))
  #current year
  ind_start <- 1
  ind_end <- records_in_week*window_weeks 
  end_week <- window_weeks
  time_no <- time_no + 1
  while (end_week <= weeks_forward)
  {
    vec_now <- vector_mean(data_now[c(ind_start:ind_end),],param_aggregation_function)
    features_now <- c(features_now,vec_now)  
    encoding_fet <- c(encoding_fet,c(1:(p-6)))
    encoding_week <- c(encoding_week,rep(1,(p-6))*time_no)
    time_no <- time_no + 1
    ind_start <- ind_start + records_in_week*step_weeks
    ind_end <- ind_end + records_in_week*step_weeks
    end_week <- ind_end/records_in_week
  }
  features_all <- rbind(features_all,features_now)
}
print(dim(features_all))
#features_all <- round(features_all,digits = 3)

sd <- apply(features_all,2,sd,na.rm = TRUE)
ind <- which(sd==0)
features_all[,ind] <- NA

#for (sk in 1:dim(features_all)[2])
#{
#  if (sum(is.na(features_all[,sk]))>0)
#  {
#    if (sum(is.na(features_all[,sk]))<16)
#    {
#      ind <- which(is.na(features_all[,sk]))
#      mn <- mean(features_all[,sk],na.rm=TRUE)
#      features_all[ind,sk] <- mn
#    }
#  }
#}

means <- apply(features_all,2,mean,na.rm=TRUE)
for (sk in 1:dim(features_all)[2])
{
  ind <- which(is.na(features_all[,sk]))
  features_all[ind,sk] <- means[sk]
}

encoding_all <- rbind(encoding_fet,encoding_week)

print(dim(features_all))
print(dim(encoding_all))


write.table(features_all,file = output_file_features_all,quote = FALSE,row.names = FALSE,col.names = FALSE,sep=',')
write.table(encoding_all,file = output_file_encoding_all,quote = FALSE,row.names = FALSE,col.names = FALSE,sep=',')
write.table(year_range,file = output_file_year_range,quote = FALSE,row.names = FALSE,col.names = FALSE,sep=',')

# table 1
names_all <- read.csv(input_file_names_all, header = TRUE, sep = '\n')

#fet_select <- c(1,2,4,8,13,15,16,17,18,21,23,27,28)
input_MS <- 'input_data/fet_MS.csv'
fet_select <- as.vector(t(read.csv(input_MS, header = FALSE)))

missing_ratio <- apply(is.na(data_all[,c(7:p)]),2,mean)
missing_ratio <- paste(as.character(ceiling(missing_ratio*100)),'\\%',sep = '')
missing_ratio[which(missing_ratio=='0\\%')] <- '-'

check <- rep(' ',length(missing_ratio))
check[fet_select] <- '\\checkmark'


table1 <- cbind(c(1:length(missing_ratio)),check,names_all,missing_ratio)

write.table(table1,file = out_table1,quote = FALSE,row.names = FALSE,col.names = FALSE,sep='&',eol='\\\\\n')
write.table(fet_select,file = out_fet_select,quote = FALSE,row.names = FALSE,col.names = FALSE)


#transformation
# plot data
for (sk in 7:dim(data_all)[2])
{
  plot_now <- paste('plots/plots_SMEAR/','fig',as.character(sk-6),'.png',sep='')
  png(plot_now,width=6,height=6,units="in",res=100)
  boxplot(data_all[,sk],main = names_all[(sk-6),])
  mn = mean(data_all[,sk],na.rm=TRUE)
  abline(h=mn,col="red")
  dev.off()
}
