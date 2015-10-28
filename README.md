Data and experiments to support a research paper **Environmental control of growth variation in a boreal Scots pine stand – a data mining approach**.
Please cite the paper if you are using the code or data in research publications.

## Detrending targets ## 


1. Detrending target variables.

		run1_detrend_target.R
	
Results are written into

	processed_data/targets_heights.csv
	processed_data/targets_widths.csv
	
	
2. Preprocessing SMEAR data.

		run2_aggregate_SMEAR.R
		
Results are written into 

	processed_data/data_aggregated_SMEAR.csv
	
Table 1 summarizing the variables is produced into

	results/table1.txt
	
3. Correlation analysis reported in Figure 2 is produced by 

		run3_correlations.R
		
Correlation tables for produciing Figure 2 are output to 

	results/corr_heights.dat	results/corr_widths.dat
	
4. Results for baseline models and traditional models

		run4_regression_greedy.R


	
	
	



