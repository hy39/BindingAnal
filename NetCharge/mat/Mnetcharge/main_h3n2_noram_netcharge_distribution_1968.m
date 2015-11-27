%main script to analyze netcharge in viral phylogenies
%1)calculate netcharge -> charge.mat
%2)merge metadata with netcharge -> merged_data.mat
%3)add gene accession to each protein accession
%4)create DNA FASTA files for BEAST input
%Oct 2, 2012
%North America 1968-2012

proj = 'h3n2_noram_1968_2012/';

%%%1)calculate charge distribution
prot_fasta = [proj 'fasta/hm_h3n2_noram_simple_1968.fas'];
charge_out = [proj 'hm_h3n2_noram_charge_1968'];
%[header charge] = main_aa_dist(prot_fasta, charge_out, 'H3N2');
disp 'Scanning amino acids';

%%%2)merge metadata with netcharge -> merged_data.mat
charge_dat = [proj 'hm_h3n2_noram_charge_1968'];
metadata = [proj 'hm_h3n2_noram_1968_ay.csv'];
merged_out = [proj 'hm_h3n2_noram_merged_data_1968'];
[v gbacc] = merge_metadata(charge_dat, metadata, merged_out);
disp 'Charged data and Date information generated.';

%%%3)add gene accession to each protein accession
merged_dat = [proj 'hm_h3n2_noram_merged_data_1968'];
accession_table_out = [proj 'hm_h3n2_noram_1968'];
%create_geneprot_table(accession_table_out, gbacc, merged_dat);
disp 'Add GenBank gene accession';

%%%4)merge with predicted binding avidity data
merged_dat = [proj 'hm_h3n2_noram_merged_data_1968'];
bding = [proj 'scr/bindingscore_h3n2_noram_noh.csv'];
merge_bding_2(bding, merged_dat);
disp 'Binding data generated';

%%%5)create DNA FASTA files for BEAST input
DNAFile = [proj 'fasta/hm_h3n2_noram_ntcds_simple_1968.fas'];
outFile = [proj 'beast/hm_h3n2_noram_beast_1968.fas'];
generate_fasta_for_beast(DNAFile, merged_dat, outFile, 1993.601, 2006.6);
disp 'BEAST input generated';