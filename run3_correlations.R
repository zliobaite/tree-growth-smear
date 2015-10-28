#2014 12 26 I.Zliobaite
#find correlation matrix

input_file_features <- 'processed_data/data_aggregated_SMEAR.csv'
input_file_encoding <- 'processed_data/encoding_aggregated_SMEAR.csv'
input_file_target_widths <- 'processed_data/targets_widths.csv'
input_file_target_heights <- 'processed_data/targets_heights.csv'
out_file_widths <- 'results/corr_widths.dat'
out_file_heights <- 'results/corr_heights.dat'
out_file_temperature <- 'results/temperature.dat'

features_all <- read.csv(input_file_features, header = FALSE, sep = ',')
encoding_all <- read.csv(input_file_encoding, header = FALSE, sep = ',')

compute_correlations <- function(features_all,encoding_all,file_target,file_out)
{
  target <- read.csv(file_target, header = FALSE, sep = '\t')
  target <- target[,2]
  
  features_all <- features_all[c(1:length(target)),]
  no_times_init <- max(encoding_all[2,])
  no_features_init <- max(encoding_all[1,])
  selection_matrix <- matrix(, nrow = (no_features_init+1), ncol = (no_times_init+1))
  
  n_fets <- dim(features_all)[2]
  
  for (sk in 1:n_fets)
  {
    CC <- cor(features_all[,sk],target)
    selection_matrix[encoding_all[1,sk],encoding_all[2,sk]] <- CC
  }
  
  selection_matrix[is.na(selection_matrix)] <- 0
  
  #writing into file
  sm <- c()
  for (sk1 in 1:dim(selection_matrix)[1])
  {
    for (sk2 in 1:dim(selection_matrix)[2])
    {
      sm <- rbind(sm,c(sk2,sk1,selection_matrix[sk1,sk2]))
    }
  }
  
  sm <- round(sm, digits = 3)
  
  write.table(sm,file = file_out,quote = FALSE,row.names = FALSE,col.names = FALSE, sep = ' ')
}


compute_correlations(features_all,encoding_all,input_file_target_widths,out_file_widths)
compute_correlations(features_all,encoding_all,input_file_target_heights,out_file_heights)