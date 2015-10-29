# 2015 02 26 I.Zliobaite
# test regressions LOOCV

input_file_features <- 'processed_data/data_aggregated_SMEAR.csv'
input_file_names <- 'input_data/names_SMEAR_short.csv'
input_file_encoding <- 'processed_data/encoding_aggregated_SMEAR.csv'
input_file_target_widths <- 'processed_data/targets_widths.csv'
input_file_target_heights <- 'processed_data/targets_heights.csv'
input_file_times <- 'input_data/week_dates.csv'

out_file_widths <- 'processed_data/greedy_R2_widths.csv'
out_file_heights <- 'processed_data/greedy_R2_heights.csv'
out_file_N2_widths <- 'processed_data/greedy_R2_N2_widths.csv'
out_file_N2_heights <- 'processed_data/greedy_R2_N2_heights.csv'
out_file_N3_widths <- 'processed_data/greedy_R2_N3_widths.csv'
out_file_N3_heights <- 'processed_data/greedy_R2_N3_heights.csv'
out_file_N4_widths <- 'processed_data/greedy_R2_N4_widths.csv'
out_file_N4_heights <- 'processed_data/greedy_R2_N4_heights.csv'
out_weeks_widths <- 'processed_data/weeks_strongest_widths.csv'
out_weeks_heights <- 'processed_data/weeks_strongest_heights.csv'
out_file_table <- 'results/table2.txt'

features_all <- read.csv(input_file_features, header = FALSE, sep = ',')
encoding_all <- read.csv(input_file_encoding, header = FALSE, sep = ',')

do_N2 <- TRUE
do_N3 <- TRUE
do_N4 <- TRUE


select_weeks <- function(features_now, encoding_now, target)
{ 
  id_features <- unique(t(encoding_now[1,]))
  id_weeks <- unique(t(encoding_now[2,]))
  
  ind_select <- c()
  cors = cor(features_now,target)
  
  cors_for_order <- c()
  for (sk in 1:length(id_features))
  {
    fet_now <- id_features[sk]
    ind_now <- which(t(encoding_now[1,])==fet_now)
    cors_now <- abs(cors[ind_now])
    max_cor <- max(cors_now,na.rm=TRUE)
    ind_sel <- which(cors_now==max_cor)
    if (length(ind_sel)>1)
    {
      ind_sel <- ind_sel[2]
      print('long index')
    }
    ind_select <- c(ind_select, ind_now[ind_sel])
  }
  return(ind_select)
}

make_regression <- function(data_training,target_training)
{
  data_training <- as.matrix(data_training)
  target_training <- as.matrix(target_training)
  beta_long = solve(t(data_training)%*%data_training)%*%(t(data_training)%*%target_training)
  return(beta_long)
}

compute_R2 <- function(true_all,predictions_all)
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

models_greedy <- function(features_all,encoding_all,file_target,do_N2,do_N3,do_N4,out_file,out_file_N2,out_file_N3,out_file_N4,out_file_week)
{
  target <- read.csv(file_target, header = FALSE, sep = '\t')
  target <- target[,2]
  yy <- length(target)
  features_all <- features_all[c(1:yy),]
  
  true_all <- c()
  predictions_N1 <- c()
  predictions_N2 <- c()
  predictions_N3 <- c()
  predictions_N4 <- c()
  week_occurence <- c()
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
    
    #extract features one week
    ind_select <- select_weeks(data_training, encoding_all, target_training)
    #print(ind_select)
    data_testing1 <- data_testing[,ind_select]
    data_training1 <- data_training[,ind_select]
    encoding1 <- encoding_all[,ind_select]
    week_occurence = rbind(week_occurence, as.matrix(encoding1[2,]))

    #make models
    pred <- c()
    for (sk2 in 1:dim(data_training1)[2])
    {
      beta_long <- make_regression(data_training1[,sk2],target_training)
      pred <- c(pred, data_testing1[,sk2]%*%beta_long + mean_target)
    }
    predictions_N1 <- rbind(predictions_N1,pred)
    
    if (do_N2)
    {
      pred <- c()
      fets_N2 <- c()
      for (sk2 in 1:(dim(data_training1)[2]-1))
      {
        for (sk3 in (sk2+1):dim(data_training1)[2])
        {
          beta_long = make_regression(data_training1[,c(sk2,sk3)],target_training)
          pred <- c(pred, as.matrix(data_testing1[,c(sk2,sk3)])%*%beta_long + mean_target)
          fets_N2 <- rbind(fets_N2,c(sk2,sk3))
        }
      }
      predictions_N2 <- rbind(predictions_N2,pred)
    }
    
    if (do_N3)
    {
      pred <- c()
      fets_N3 <- c()
      #print(dim(data_training1))
      for (sk2 in 1:(dim(data_training1)[2]-2))
      {
        for (sk3 in (sk2+1):(dim(data_training1)[2]-1))
        {
          for (sk4 in (sk3+1):dim(data_training1)[2])
          {
            beta_long = make_regression(data_training1[,c(sk2,sk3,sk4)],target_training)
            pred <- c(pred, as.matrix(data_testing1[,c(sk2,sk3,sk4)])%*%beta_long + mean_target)
            fets_N3 <- rbind(fets_N3,c(sk2,sk3,sk4))  
          }
        }
      }
      predictions_N3 <- rbind(predictions_N3,pred)
    }
    
    if (do_N4)
    {
      pred <- c()
      fets_N4 <- c()
      for (sk2 in 1:(dim(data_training1)[2]-3))
      {
        for (sk3 in (sk2+1):(dim(data_training1)[2]-2))
        {
          for (sk4 in (sk3+1):(dim(data_training1)[2]-1))
          {
            for (sk5 in (sk4+1):dim(data_training1)[2])
            {
              beta_long = make_regression(data_training1[,c(sk2,sk3,sk4,sk5)],target_training)
              pred <- c(pred, as.matrix(data_testing1[,c(sk2,sk3,sk4,sk5)])%*%beta_long + mean_target)
              fets_N4 <- rbind(fets_N4,c(sk2,sk3,sk4,sk5))  
            }
          }
        }
      }
      predictions_N4 <- rbind(predictions_N4,pred)
    }
  }
  
  #results
  R2 <- compute_R2(true_all,predictions_N1)
  results <- cbind(c(1:length(R2)),R2)
  if (do_N2)
  {
    R2 <- compute_R2(true_all,predictions_N2)
    results_N2 <- cbind(fets_N2,R2)
  }
  if (do_N3)
  {
    R2 <- compute_R2(true_all,predictions_N3)
    results_N3 <- cbind(fets_N3,R2)
  }
  if (do_N4)
  {
    R2 <- compute_R2(true_all,predictions_N4)
    results_N4 <- cbind(fets_N4,R2)
  }
  
  #training fit
  mean_now <- apply(features_all,2,mean,na.rm=TRUE)
  std_now <- apply(features_all,2,sd,na.rm=TRUE)
  std_now[std_now==0] <- 1
  features_new <- features_all - matrix(1,yy,1)%*%t(as.matrix(mean_now))
  features_new <- features_new/matrix(1,yy,1)%*%t(as.matrix(std_now))
  
  #replace missing values by mean (0)
  features_new[is.na(features_new)] <- 0
  
  mean_target <- mean(target)
  target <- target - mean_target
  #extract features one week
  ind_select <- select_weeks(features_new, encoding_all, target)
  features_new <- features_new[,ind_select]
  encoding_new <- encoding_all[,ind_select]
  week_occurence = rbind(week_occurence, as.matrix(encoding_new[2,]))
  
  #make models
  predictions_fit <- c()
  for (sk2 in 1:dim(features_new)[2])
  {
    beta_long <- make_regression(features_new[,sk2],target)
    pred <- features_new[,sk2]%*%beta_long + mean_target
    predictions_fit <- cbind(predictions_fit, pred)
  }
  R2 <- compute_R2(true_all,predictions_fit)
  results <- cbind(results,R2)
  colnames(results) <- c('id','R2_N1_cv','R2_N1_fit')
  
  if (do_N2)
  {
    predictions_fit <- c()
    for (sk2 in 1:(dim(features_new)[2]-1))
    {
      for (sk3 in (sk2+1):dim(features_new)[2])
      {
        beta_long <- make_regression(features_new[,c(sk2,sk3)],target)
        pred <- as.matrix(features_new[,c(sk2,sk3)])%*%beta_long + mean_target
        predictions_fit <- cbind(predictions_fit, pred)  
      }
    }
    R2 <- compute_R2(true_all,predictions_fit)
    results_N2 <- cbind(results_N2,R2)
    colnames(results_N2) <- c('id1','id2','R2_N2_cv','R2_N2_fit')
  }
  
  if (do_N3)
  {
    predictions_fit <- c()
    for (sk2 in 1:(dim(features_new)[2]-2))
    {
      for (sk3 in (sk2+1):(dim(features_new)[2]-1))
      {
        for (sk4 in (sk3+1):dim(features_new)[2])
        {
          beta_long <- make_regression(features_new[,c(sk2,sk3,sk4)],target)
          pred <- as.matrix(features_new[,c(sk2,sk3,sk4)])%*%beta_long + mean_target
          predictions_fit <- cbind(predictions_fit, pred)  
        }
      }
    }
    R2 <- compute_R2(true_all,predictions_fit)
    results_N3 <- cbind(results_N3,R2)
    colnames(results_N3) <- c('id1','id2','id3','R2_N3_cv','R2_N3_fit')
  }
  
  if (do_N4)
  {
    predictions_fit <- c()
    for (sk2 in 1:(dim(features_new)[2]-3))
    {
      for (sk3 in (sk2+1):(dim(features_new)[2]-2))
      {
        for (sk4 in (sk3+1):(dim(features_new)[2]-1))
        {
          for (sk5 in (sk4+1):dim(features_new)[2])
          {
            beta_long <- make_regression(features_new[,c(sk2,sk3,sk4,sk5)],target)
            pred <- as.matrix(features_new[,c(sk2,sk3,sk4,sk5)])%*%beta_long + mean_target
            predictions_fit <- cbind(predictions_fit, pred)  
          }
        }
      }
    }
    R2 <- compute_R2(true_all,predictions_fit)
    results_N4 <- cbind(results_N4,R2)
    colnames(results_N4) <- c('id1','id2','id3','id4','R2_N4_cv','R2_N4_fit')
  }
  
  write.table(round(results,digits = 3),file = out_file,quote = FALSE,row.names = FALSE,col.names = TRUE, sep = ',')
  write.table(week_occurence,file = out_file_week,quote = FALSE,row.names = FALSE,col.names = FALSE, sep = ',')
  write.table(round(results_N2,digits = 3),file = out_file_N2,quote = FALSE,row.names = FALSE,col.names = TRUE, sep = ',')
  write.table(round(results_N3,digits = 3),file = out_file_N3,quote = FALSE,row.names = FALSE,col.names = TRUE, sep = ',')
  write.table(round(results_N4,digits = 3),file = out_file_N4,quote = FALSE,row.names = FALSE,col.names = TRUE, sep = ',')
}


models_greedy(features_all,encoding_all,input_file_target_widths,do_N2,do_N3,do_N4,out_file_widths,out_file_N2_widths,out_file_N3_widths,out_file_N4_widths,out_weeks_widths)
models_greedy(features_all,encoding_all,input_file_target_heights,do_N2,do_N3,do_N4,out_file_heights,out_file_N2_heights,out_file_N3_heights,out_file_N4_heights,out_weeks_heights)


#write results, table 4
n_width <- 16
n_height <- 15  
n_fets <- 31

mode_fun <- function(x) 
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

all_modes_fun <- function(data1)
{
  md <- c()
  for (sk in 1:dim(data1)[2])
  {
    md <- c(md,mode_fun(data1[,sk]))
  }
  return(md)
}

stability_fun <- function(data1,modes)
{
  st <- c()
  n <- dim(data1)[1]
  for (sk in 1:length(modes))
  {
    st <- c(st,sum(data1[,sk]==modes[sk])/n)
  }
  return(round(st,digits=1))
}

times_fun <- function(times,modes)
{
  n <- length(modes)
  tm <- c()
  for (sk in 1:n)
  {
    tm <- c(tm, as.vector(times[modes[sk],]))
  }
  return(tm)
}

#read names
fet_names <- read.table(input_file_names, header = TRUE, '\t')

#read time
week_times <- read.table(input_file_times, header = TRUE, '\t')

#read weeks
weeks_width <- read.csv(out_weeks_widths, header = FALSE, sep = ',')
modes_width <- all_modes_fun(weeks_width[c(1:n_width),])
stab_width <- stability_fun(weeks_width[c(1:n_width),],modes_width)
times_width <- times_fun(week_times,modes_width)

weeks_height <- read.csv(out_weeks_heights, header = FALSE, sep = ',')
modes_height <- all_modes_fun(weeks_height[c(1:n_height),])
stab_height <- stability_fun(weeks_height[c(1:n_width),],modes_height)
times_height <- times_fun(week_times,modes_height)

#read R2
R2_width <- read.csv(out_file_widths, header = TRUE, sep = ',')
R2_width <- round(R2_width,digits = 2)

R2_height <- read.csv(out_file_heights, header = TRUE, sep = ',')
R2_height <- round(R2_height,digits = 2)

#as.vector(t(weeks_width[n_width+1,])) #main model weeks (the same)
table_all <- cbind(c(1:n_fets),fet_names,times_width,stab_width,R2_width[,3],R2_width[,2],times_height,stab_height,R2_height[,3],R2_height[,2])
colnames(table_all) <- c('No','Feature','Time','Cons.','R2 fit','R2 test','Time','Cons.','R2 fit','R2 test')

write.table(table_all,file = out_file_table,quote = FALSE,row.names = FALSE,sep=' & ',eol = '\\\\ \n')
