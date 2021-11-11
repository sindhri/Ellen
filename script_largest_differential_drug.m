%20210617, specify target and to_exclude_LL and run the pipeline
%find largest differential effect
%1. calculate z-score differentially against WT+DMSO and HOM+DMSO
%2. averege the z-scores across fish
%3. calcualte euclidean distance
%4. pair HOM/HET and WT for the same drug
%5. calculate the averaged differential effect for each drug
addpath('helpers/');
largest_differential_drug('HOM', 0);
largest_differential_drug('HOM', 1);

largest_differential_drug('HET', 0);
largest_differential_drug('HET', 1);

