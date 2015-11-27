%main script to analyze netcharge in viral phylogenies
%1)calculate netcharge -> charge.mat
%2)merge metadata with netcharge -> merged_data.mat
%3)add gene accession to each protein accession
%4)create DNA FASTA files for BEAST input
%Sep 26, 2012
%New Zealand 1999-2005

proj = 'h3n2_sea/';
prot_fasta = [proj 'hm_h3n2_sea_any_simple.fas'];
charge_out = [proj 'hm_h3n2_sea_charge'];
[header charge] = main_aa_dist(prot_fasta, charge_out, 'H3N2');

charge_dat = [proj 'hm_h3n2_sea_charge'];
metadata = [proj 'hm_h3n2_sea_any_ay.csv'];
merged_out = [proj 'hm_h3n2_sea_merged_data'];
[v gbacc] = merge_metadata(charge_dat, metadata, merged_out);

merged_dat = [proj 'hm_h3n2_sea_merged_data'];
accession_table_out = [proj 'hm_h3n2_sea_any'];
create_geneprot_table(accession_table_out, gbacc, merged_dat);
disp 'Charged data generated';

DNAFile = [proj 'hm_h3n2_sea_ntcds_simple.fas'];
outFile = [proj 'hm_h3n2_sea_dna_beast_1993_2006.fas'];
generate_fasta_for_beast(DNAFile, merged_dat, outFile, 1993.601, 2006.6);
disp 'BEAST input generated';