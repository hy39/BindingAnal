%main script to analyze netcharge in viral phylogenies
%2)calculate netcharge -> charge.mat
%3)merge metadata with netcharge -> merged_data.mat
%4)add gene accession to each protein accession
%5)create DNA FASTA files for BEAST input
%Sep 26, 2012

proj = 'h1n1_world_1918_2007/';

%step1 exclude outliers
metadata = [proj 'hm_h1n1_whole_ay.csv'];
metadata_exo = [proj 'hm_h1n1_flu_world_ay_exo.csv'];
exc_dat = [proj 'excluded_strains_2.txt'];
%label_excluded_data(metadata, metadata_exo, exc_dat);
%save into metadata_exo

%step2calculate netcharge
%prot_fasta = [proj 'hm_h1n1_flu_whole_simple.fas'];
prot_fasta = [proj 'fasta/hm_h1n1_world_aligned.fas'];
charge_out = [proj 'hm_h1n1_whole_charge'];
%[header charge] = main_aa_dist(prot_fasta, charge_out, 'H1N1');

%step3 merge data, netcharge and ages information
charge_dat = [proj 'hm_h1n1_whole_charge'];
metadata = [proj 'hm_h1n1_flu_world_ay_exo.csv']; %use metadata_exo
%%metadata = [proj 'hm_h1n1_whole_ay.csv'];
merged_out = [proj 'hm_h1n1_whole_merged_data'];
%[v gbacc] = merge_metadata(charge_dat, metadata, merged_out);

%step4 merge with predicted binding avidity data
merged_dat = [proj 'hm_h1n1_whole_merged_data'];
bding = [proj 'scr/bindingscore_h1n1_full_noh.csv'];
%merge_bding_2(bding, merged_dat);
disp 'Binding data generated';

%step5 add gene accession
merged_dat = [proj 'hm_h1n1_whole_merged_data'];
%accession_table = [proj 'hm_h1n1_flu_whole'];
accession_table = [proj 'hm_h1n1'];
%create_geneprot_table(accession_table, gbacc, merged_dat);
%disp 'add both gene and protein accession';

%step6 generate FASTA file
DNAFile = [proj 'hm_h1n1_noram_ntcds_simple1.fas'];
outFile = [proj 'hm_h1n1_noram_age_dna_beast_1995_2008.fas'];
generate_fasta_for_beast(DNAFile, merged_dat, outFile, 1995, 2008.999);
%disp 'BEAST input generated';