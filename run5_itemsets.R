# 2015 04 05 I.Zliobaite

input_file_widths <- 'processed_data/greedy_R2_widths.csv'
input_file_heights <- 'processed_data/greedy_R2_heights.csv'
input_file_N2_widths <- 'processed_data/greedy_R2_N2_widths.csv'
input_file_N2_heights <- 'processed_data/greedy_R2_N2_heights.csv'
input_file_N3_widths <- 'processed_data/greedy_R2_N3_widths.csv'
input_file_N3_heights <- 'processed_data/greedy_R2_N3_heights.csv'
input_file_N4_widths <- 'processed_data/greedy_R2_N4_widths.csv'
input_file_N4_heights <- 'processed_data/greedy_R2_N4_heights.csv'
input_file_names <- 'input_data/names_SMEAR_short.csv'

out_file_table_a <- 'results/table3a.txt'
out_file_table_b <- 'results/table3b.txt'

param_no_select <- 15
param_in_thresh <- 0.15

update_indicator <- function(combinations,indicator,ind)
{
    comb_now <- combinations[ind,]
    comb_now <- comb_now[comb_now>0]
    for (sk2 in 1:length(indicator))
    {
      comp_now <- combinations[sk2,]
      comp_now <- comp_now[comp_now>0]
      inter <- intersect(comp_now,comb_now)
      if (length(comb_now)==1)
      {
        if (length(inter)>0)
        {
          indicator[sk2] <- 0
        }
      }else
      {
        if (length(inter)>1)
        {
          indicator[sk2] <- 0
        }
      }
    }
  return(indicator)
}


select_subsets <- function(no_select,file_N1,file_N2,file_N3,file_N4)
{
  data_N1 <- read.csv(file_N1,header = TRUE, sep = ',')
  data_N1 <- cbind(data_N1[,1],matrix(0,dim(data_N1)[1],3),data_N1[,c(2,3)])
  data_N2 <- read.csv(file_N2,header = TRUE, sep = ',')
  data_N2 <- cbind(data_N2[,c(1,2)],matrix(0,dim(data_N2)[1],2),data_N2[,c(3,4)])
  data_N3 <- read.csv(file_N3,header = TRUE, sep = ',')
  data_N3 <- cbind(data_N3[,c(1,2,3)],matrix(0,dim(data_N3)[1],1),data_N3[,c(4,5)])
  data_N4 <- read.csv(file_N4,header = TRUE, sep = ',')
  data_all <- rbind(as.matrix(data_N1),as.matrix(data_N2),as.matrix(data_N3),as.matrix(data_N4))
  colnames(data_all) <- c('id1','id2','id3','id4','R2cv','R2fit')
  indicator <- matrix(1,dim(data_all)[1],1)
  selected <- c()
  ind_open <- which(indicator==1)
  for (sk in 1:no_select)
  {
    ind <- which(data_all[ind_open,5]==max(data_all[ind_open,5]))
    ind <- ind_open[ind]
    if(length(ind)>1)
    {
      ln <-0
      for (sk2 in 1:length(ind))
      {
        ln <- c(ln,length(which(data_all[ind[sk],c(1:4)]>0)))
      }
      ind_sel <- which(ln==min(ln))
      if (length(ind_sel)>1){ind_sel<- ind_sel[1]}
      ind <- ind[ind_sel]
    }
    data_now <- data_all[ind,c(1:4)]
    if (length(which(data_now>0))>1)
    {
      if (data_all[ind,'R2cv']>=param_in_thresh)
      {
        selected <- rbind(selected,data_all[ind,])    
      }
    }
    indicator <- update_indicator(data_all[,c(1:4)],indicator,ind)
    ind_open <- which(indicator==1)
  }
  return(selected)
}

ids_to_names <- function(data_selected,names)
{
  table_all <- c()
  for (sk in 1:dim(data_selected)[1])
  {
    name_row <- c()
    for (sk2 in 1:4)
    {
      if (data_selected[sk,sk2] == 0)
      {
        name_now = ' '
      }else
      {
        name_now = as.matrix(names[data_selected[sk,sk2],])
      }
      name_row <- cbind(name_row,name_now)
    }
    table_all <- rbind(table_all,c(name_row,round(data_selected[sk,c(5,6)],digits=2)))
  }
  return(table_all)
}

data_selected_widths <- select_subsets(param_no_select,input_file_widths,input_file_N2_widths,input_file_N3_widths,input_file_N4_widths)
data_selected_heights <- select_subsets(param_no_select,input_file_heights,input_file_N2_heights,input_file_N3_heights,input_file_N4_heights)

#write results
names_all <- read.csv(input_file_names, header = TRUE, sep = '\t')
table1 <- ids_to_names(data_selected_widths,names_all)
table2 <- ids_to_names(data_selected_heights,names_all)

write.table(table1,file = out_file_table_a,quote = FALSE,row.names = FALSE,col.names = FALSE,sep=' & ',eol = '\\\\ \n')
write.table(table2,file = out_file_table_b,quote = FALSE,row.names = FALSE,col.names = FALSE,sep=' & ',eol = '\\\\ \n')