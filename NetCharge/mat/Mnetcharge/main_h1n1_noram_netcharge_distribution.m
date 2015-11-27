%main script to analyze netcharge in viral phylogenies
%2)calculate netcharge -> charge.mat
%3)merge metadata with netcharge -> merged_data.mat
%4)add gene accession to each protein accession
%5)create DNA FASTA files for BEAST input
%Sep 26, 2012

p = path
path(p,'lib/');

proj = 'h1n1_noram_1995/';

%step1 exclude outliers
metadata = [proj 'hm_h1n1_flu_noram_ay_new.csv'];
metadata_exo = [proj 'hm_h1n1_flu_noram_ay_exo.csv']; %metadata output
exc_dat = [proj 'excluded_strains_2.txt'];
label_excluded_data(metadata, metadata_exo, exc_dat);
disp 'exclude outliers';

%step2 calculate netcharge
prot_fasta = [proj 'fasta/hm_h1n1_flu_noram_any_simple_new.fas'];
charge_out = [proj 'hm_h1n1_noram_charge'];
[header charge] = main_aa_dist(prot_fasta, charge_out, 'H1N1');
disp 'calculate netcharge';

%step3 merge data, netcharge and ages information
charge_dat = [proj 'hm_h1n1_noram_charge'];
metadata = [proj 'hm_h1n1_flu_noram_ay_exo.csv'];
merged_out = [proj 'hm_h1n1_noram_merged_data'];
[v gbacc] = merge_metadata(charge_dat, metadata, merged_out);

%step4 merge with predicted binding avidity data
merged_dat = [proj 'hm_h1n1_noram_merged_data'];
bding = [proj 'scr/bindingscore_h1n1_noram_1995_2008.csv'];
merge_bding_2(bding, merged_dat);
disp 'Binding data generated';

%step5 add gene accession
merged_dat = [proj 'hm_h1n1_noram_merged_data'];
accession_table = [proj 'hm_h1n1_flu_noram'];
create_geneprot_table(accession_table, gbacc, merged_dat);
disp 'Charged data generated';

%step6 generate FASTA file
DNAFile = [proj 'fasta/hm_h1n1_noram_ntcds_simple1_new.fas'];
outFile = [proj 'hm_h1n1_noram_age_dna_beast_1995_2006.fas'];
%generate_fasta_for_beast(DNAFile, merged_dat, outFile, 1995, 2008.999);
generate_fasta_for_beast_full(DNAFile, merged_dat, outFile, 1995, 2006.5);
disp 'BEAST input generated';