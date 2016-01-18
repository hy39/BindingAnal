%main script to analyze netcharge in viral phylogenies
%1)calculate netcharge -> charge.mat
%2)merge metadata with netcharge -> merged_data.mat
%3)add gene accession to each protein accession
%4)create DNA FASTA files for BEAST input
%Sep 26, 2012

prot_fasta = 'h1n1/hm_h1n1_flu_noram_any_simple.fas';
charge_out = 'h1n1/hm_h1n1_noram_charge';
[header charge] = main_aa_dist(prot_fasta, charge_out);

charge_dat = 'h1n1/hm_h1n1_noram_charge';
metadata = 'h1n1/hm_h1n1_flu_noram_ay.csv';
merged_out = 'h1n1/hm_h1n1_noram_merged_data';
[v gbacc] = merge_metadata(charge_dat, metadata, merged_out);

merged_dat = 'h1n1/hm_h1n1_noram_merged_data';
accession_table = 'h1n1/hm_h1n1_flu_noram';
create_geneprot_table(accession_table, gbacc, merged_dat);
disp 'Charged data generated';

DNAFile = 'h1n1/hm_h1n1_noram_ntcds_simple1.fas';
outFile = 'h1n1/hm_h1n1_noram_age_dna_beast_2000_2008.fas';
format_fasta_header(DNAFile, merged_dat, outFile);
disp 'BEAST input generated';