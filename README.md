Data and experiments to support a research paper **Environmental control of growth variation in a boreal Scots pine stand – a data mining approach**.
Please cite the paper if you are using the code or data in research publications.

Paper: *Liisa Kulmala, Indre Zliobaite, Eero Nikinmaa, Pekka Nojd, Pasi Kolari, K. Kabiri Koupaei, Jaakko Hollmen, Harri Makinen (2015). Environmental control of growth variation in a boreal Scots pine stand – a data mining approach. Under review.*

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
		
Correlation tables for producing Figure 2 are output to 

	results/corr_heights.dat	results/corr_widths.dat
		
4. Results for baseline models and traditional models

		run4_regression_greedy.R
		
Produces Table 2 with accuracies for single variable models

	results/table2.txt

5. Results for greedy approach 

		run5_itemsets.R
		
Produces Table 3		
	
	results/table3a.txt	
	results/table3b.txt	
	
6. LARS regression

		run6_regression_LARS.R
		
Produces accuracy records for making Figure 3,

	results/out_R2_LARS_heights.dat
	results/out_R2_LARS_widths.dat
	

and Figure 4

	results/cofs_LARS_BB_heights.dat
	results/cofs_LARS_BB_widths.dat
	results/cofs_LARS_MS_heights.dat
	results/cofs_LARS_MS_widths.dat
	results/cofs_LARS_TR_heights.dat
	results/cofs_LARS_TR_widths.dat



