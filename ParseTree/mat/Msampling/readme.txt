% Collect influenza sequences whole set. 1968 – 2013/07/01
fasta2mat('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds.fas'); 
generate_random_years('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds.mat',305,1968,2013)
generate_random_tree('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds.mat', 'ys_305_2.mat')


% Generate standard set
removeOutlier('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_outlier_sets.csv', 'dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds.fas');
fasta2mat('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds_std.fas'); 
generate_random_years('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds_std.mat',300,1968,2013)
generate_random_tree('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds_std.mat', 'ys_300_3.mat')

###
Notes:
The records before 1990 (or maybe 1991) lack of month and date information.
Assign the date Jul 1 as month and day information.