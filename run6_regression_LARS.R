#2015 03 10 I.Zliobaite
#LARS regression on everything

library(lars)

input_MS <- 'input_data/fet_MS.csv'
input_TR <- 'input_data/fet_TR.csv'

input_file_features <- 'processed_data/data_aggregated_SMEAR.csv'
input_file_encoding <- 'processed_data/encoding_aggregated_SMEAR.csv'
input_file_target_widths <- 'processed_data/targets_widths.csv'
input_file_target_heights <- 'processed_data/targets_heights.csv'

out_file_name_BB_wd = 'results/cofs_LARS_BB_widths.dat'
out_file_name_MS_wd = 'results/cofs_LARS_MS_widths.dat'
out_file_name_TR_wd = 'results/cofs_LARS_TR_widths.dat'
out_file_R2_wd = 'results/out_R2_LARS_widths.dat'

out_file_name_BB_ht = 'results/cofs_LARS_BB_heights.dat'
out_file_name_MS_ht = 'results/cofs_LARS_MS_heights.dat'
out_file_name_TR_ht = 'results/cofs_LARS_TR_heights.dat'
out_file_R2_ht = 'results/out_R2_LARS_heights.dat'

out_file_tmp = 'results/out_mean_tmp.dat'

param_steps <- 10

features_all <- read.csv(input_file_features, header = FALSE, sep = ',')
encoding_all <- read.csv(input_file_encoding, header = FALSE, sep = ',')

fet_MS <- read.csv(input_MS, header = FALSE)
fet_TR <- read.csv(input_TR, header = FALSE)


expand_ind <- function(encoding_all,ind)
{
  ind_all <- c()
  for (sk in 1:dim(ind)[1])
  {
    ind_now <- which(encoding_all[1,]==ind[sk])
    ind_all <- c(ind_all,ind_now)
  }
  return(ind_all)
}

compute_errors <- function(true_all,predictions_all)
{
  p <- dim(predictions_all)[2]  
  R2 <- c()
  for (sk in 1:p)
  {
    R2_now <- 1 - t(true_all - predictions_all[,sk])%*%(true_all - predictions_all[,sk])/(t(true_all - mean(true_all))%*%(true_all - mean(true_all)))
    R2 <- rbind(R2,R2_now)
  }
  return(R2)
}

compute_cof_stability <- function(beta)
{
  mn <- apply(beta,2,mean)
  ln <- length(beta)
  arr <- mn*ln/100  
}

make_plot_matrix <- function(arr,encoding_all)
{
  n_cols <- length(unique(t(encoding_all[2,])))  
  n_rows <- length(arr)/n_cols
  enco <- t(encoding_all[1,])
  fets <- unique(enco)
  if (max(fets)>length(fets))
  {
    for (sk in 1:length(fets))
    {
      ind <- which(encoding_all[1,] == fets[sk])
      enco[ind] <- sk
    }
  }
  
  cofs <- matrix(0,n_rows+1,n_cols+1)
  
  for (sk in 1:length(arr))
  {
    cofs[enco[sk],encoding_all[2,sk]] <- arr[sk]
  }
  
  #for writing into file
  sm <- c()
  for (sk1 in 1:(n_rows+1))
  {
    for (sk2 in 1:(n_cols+1))
    {
      sm <- rbind(sm, c(sk2,sk1,cofs[sk1,sk2]))  
    }    
  }
  return(sm)
}


compute_LARS <- function(features_all,encoding_all,file_target,fet_MS,fet_TR,out_file_name_BB,out_file_name_MS,out_file_name_TR,out_file_R2)
{
  target <- read.csv(file_target, header = FALSE, sep = '\t')
  target <- target[,2]
  
  features_all <- features_all[c(1:length(target)),]
  yy <- dim(features_all)[1]
  
  fet_MS = expand_ind(encoding_all,t(fet_MS))
  fet_TR = expand_ind(encoding_all,t(fet_TR))

  true_all <- c()
  predictions_BB <- c()
  predictions_MS <- c()
  predictions_TR <- c()
  predictions_nai <- c()
  predictions_BB_fit <- c()
  predictions_MS_fit <- c()
  predictions_TR_fit <- c()
  Ball_BB <- c()
  Ball_MS <- c()
  Ball_TR <- c()
  for (ind_testing in 1:yy)
  {
    #print(ind_testing)
    #training and testing indices
    ind_training <- setdiff(c(1:yy),ind_testing)
  
    #training and testing data
    data_training <- features_all[ind_training,]
    data_testing <- features_all[ind_testing,]
    target_training <- target[ind_training]
    target_testing <- target[ind_testing]
    mean_target <- mean(target_training)
    target_training <- target_training - mean_target
    true_all <- c(true_all, target_testing)
  
    #normalize data
    mean_now <- apply(data_training,2,mean,na.rm=TRUE)
    std_now <- apply(data_training,2,sd,na.rm=TRUE)
    std_now[std_now==0] <- 1
    data_training <- data_training - matrix(1,length(ind_training),1)%*%t(as.matrix(mean_now))
    data_training <- data_training/matrix(1,length(ind_training),1)%*%t(as.matrix(std_now))
    data_testing <- data_testing - mean_now
    data_testing <- data_testing/std_now
  
    #replace missing values by mean (0)
    data_training[is.na(data_training)] <- 0
    data_testing[is.na(data_testing)] <- 0
  
    # LARS regression
    tp <- 'lar'
    nml <- FALSE
    fit_lars <- lars(as.matrix(data_training),as.matrix(target_training), normalize = nml, max.steps = param_steps, use.Gram=FALSE, type = tp, intercept = FALSE)
    model_now <- coef(fit_lars, mode = 'step')
    model_now <- model_now[11,]
    Ball_BB <- rbind(Ball_BB,model_now)
    intercepts_lars <- predict.lars(fit_lars, data_training*0, type="fit")
    pred <- predict.lars(fit_lars, data_testing, type="fit")
    predictions_BB <- rbind(predictions_BB,(pred$fit + mean_target))
  
    
    fit_lars <- lars(as.matrix(data_training[,fet_MS]),as.matrix(target_training), normalize = nml, max.steps = param_steps, use.Gram=FALSE, type = tp, intercept = FALSE)
    model_now <- coef(fit_lars, mode = 'step')
    model_now <- model_now[11,]
    Ball_MS <- rbind(Ball_MS,model_now)
    intercepts_lars <- predict.lars(fit_lars, data_training[,fet_MS]*0, type="fit")
    pred <- predict.lars(fit_lars, data_testing[,fet_MS], type="fit")
    predictions_MS <- rbind(predictions_MS,(pred$fit + mean_target))
  
    fit_lars <- lars(as.matrix(data_training[,fet_TR]),as.matrix(target_training), normalize = nml, max.steps = param_steps, use.Gram=FALSE, type = tp,intercept = FALSE)
    model_now <- coef(fit_lars, mode = 'step')
    model_now <- model_now[11,]
    Ball_TR <- rbind(Ball_TR,model_now)
    intercepts_lars <- predict.lars(fit_lars, data_training[,fet_TR]*0, type="fit")
    pred <- predict.lars(fit_lars, data_testing[,fet_TR], type="fit")
    predictions_TR <- rbind(predictions_TR,(pred$fit + mean_target))

    predictions_nai <- rbind(predictions_nai,(mean(target_training)+mean_target))
  }

  #training fit
  mean_now = apply(features_all,2,mean,na.rm=TRUE)
  std_now = apply(features_all,2,sd,na.rm=TRUE)
  std_now[std_now==0] <- 1
  features_new <- features_all - matrix(1,yy,1)%*%t(as.matrix(mean_now))
  features_new <- features_new/matrix(1,yy,1)%*%t(as.matrix(std_now))

  mean_target <- mean(target)
  target <- target - mean_target

  #replace missing values by mean (0)
  features_new[is.na(features_new)] <- 0

  fit_lars <- lars(as.matrix(features_new),as.matrix(target), normalize = nml, max.steps = param_steps, use.Gram=FALSE, type = tp, intercept = FALSE)
  intercepts_lars <- predict.lars(fit_lars, features_new*0, type="fit")
  pred <- predict.lars(fit_lars, features_new, type="fit")
  predictions_BB_fit <- rbind(predictions_BB_fit,(pred$fit + mean_target))

  fit_lars <- lars(as.matrix(features_new[,fet_MS]),as.matrix(target), normalize = nml, max.steps = param_steps, use.Gram=FALSE, type = tp, intercept = FALSE)
  intercepts_lars <- predict.lars(fit_lars, features_new[,fet_MS]*0, type="fit")
  pred <- predict.lars(fit_lars, features_new[,fet_MS], type="fit")
  predictions_MS_fit <- rbind(predictions_MS_fit,(pred$fit + mean_target))

  fit_lars <- lars(as.matrix(features_new[,fet_TR]),as.matrix(target), normalize = nml, max.steps = param_steps, use.Gram=FALSE, type = tp, intercept = FALSE)
  intercepts_lars <- predict.lars(fit_lars, features_new[,fet_TR]*0, type="fit")
  pred <- predict.lars(fit_lars, features_new[,fet_TR], type="fit")
  predictions_TR_fit <- rbind(predictions_TR_fit,(pred$fit + mean_target))

  R2_BB <- compute_errors(true_all,predictions_BB)
  R2_MS <- compute_errors(true_all,predictions_MS)
  R2_TR <- compute_errors(true_all,predictions_TR)
  R2_nai <- compute_errors(true_all,predictions_nai)

  R2_BB_fit <- compute_errors(true_all,predictions_BB_fit)
  R2_MS_fit <- compute_errors(true_all,predictions_MS_fit)
  R2_TR_fit <- compute_errors(true_all,predictions_TR_fit)

  results <- cbind(c(0:(length(R2_BB)-1)),R2_BB_fit,R2_BB,R2_MS_fit,R2_MS,R2_TR_fit,R2_TR,matrix(1,length(R2_TR),1)%*%R2_nai)
  colnames(results) <- c('nn','R2BBfit','R2BB','R2MSfit','R2MS','R2TRfit','R2TR','base')
  write.table(round(results,digits = 3), file = out_file_R2, row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)

  Ball_BB = compute_cof_stability(Ball_BB) 
  sm = make_plot_matrix(Ball_BB,encoding_all)
  write.table(sm, file = out_file_name_BB, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

  #print('BB done')
  Ball_MS = compute_cof_stability(Ball_MS) 
  sm = make_plot_matrix(Ball_MS,encoding_all[,fet_MS])
  write.table(sm, file = out_file_name_MS, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

  #print('MS done')
  Ball_TR = compute_cof_stability(Ball_TR) 
  sm = make_plot_matrix(Ball_TR,encoding_all[,fet_TR])
  write.table(sm, file = out_file_name_TR, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
  
}

compute_LARS(features_all,encoding_all,input_file_target_widths,fet_MS,fet_TR,out_file_name_BB_wd,out_file_name_MS_wd,out_file_name_TR_wd,out_file_R2_wd)
compute_LARS(features_all,encoding_all,input_file_target_heights,fet_MS,fet_TR,out_file_name_BB_ht,out_file_name_MS_ht,out_file_name_TR_ht,out_file_R2_ht)

#save mean temperature
un_weeks <- as.vector(unique(t(encoding_all[2,])))
mean_tmp <- c()
for (sk in 1:length(un_weeks))
{
  ind <- intersect(which(encoding_all[1,]==2),which(encoding_all[2,]==un_weeks[sk]))
  mean_tmp <- c(mean_tmp,mean(features_all[,ind]))
}
mean_tmp <- cbind(un_weeks,round(mean_tmp,digits = 2))
colnames(mean_tmp) <- c('wk','tmp')
write.table(mean_tmp, file = out_file_tmp, row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)