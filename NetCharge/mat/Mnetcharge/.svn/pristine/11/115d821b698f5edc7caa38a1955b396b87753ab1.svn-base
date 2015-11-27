%main script to analyze netcharge in viral phylogenies
%1)calculate netcharge -> charge.mat
%2)merge metadata with netcharge -> merged_data.mat
%3)add gene accession to each protein accession
%4)create DNA FASTA files for BEAST input
%v1, Oct 2, 2012
%v2, Feb 4, 2013
%New York 1993.6-2006.6

%p = path
%path(p,'../');
p = path
path(p,'lib/');

proj = 'h3n2_ny/';

%%%1)calculate charge distribution
prot_fasta = [proj 'fasta/hm_h3n2_ny_any_simple.fas'];
charge_out = [proj 'hm_h3n2_ny_charge'];
%[header charge] = main_aa_dist(prot_fasta, charge_out, 'H3N2');
disp 'Scanning amino acids';

%%%2)merge metadata with netcharge -> merged_data.mat
charge_dat = [proj 'hm_h3n2_ny_charge'];
metadata = [proj 'hm_h3n2_ny_ay.csv'];
redundant_set = [proj 'hm_h3n2_ny_redundant_sets.csv'];
metadata_nr_out = [proj 'hm_h3n2_ny_ay_nr.csv'];
merged_out = [proj 'hm_h3n2_ny_merged_data'];
%% For all data, might contains redundant
[v gbacc] = merge_metadata(charge_dat, metadata, merged_out);
%% For non-redundant set
%%[v gbacc] = merge_metadata_nr(charge_dat, metadata, redundant_set, merged_out);
disp 'Charged data generated';

%%%3)add gene accession to each protein accession
merged_dat = [proj 'hm_h3n2_ny_merged_data'];
accession_table_out = [proj 'hm_h3n2_ny_any'];
create_geneprot_table(accession_table_out, gbacc, merged_dat); %add gene accession
%%convert_age2k(merged_dat);
disp 'Add GenBank gene accession';

%%%4)merge with predicted binding avidity data
merged_dat = [proj 'hm_h3n2_ny_merged_data'];
proj1 = 'h3n2_noram_1968_2012/';
bding = [proj1 'scr/bindingscore_h3n2_noram_noh.csv'];
merge_bding_2(bding, merged_dat);
disp 'Binding data generated';

%%%5)create DNA FASTA files for BEAST input
merged_dat = [proj 'hm_h3n2_ny_merged_data'];
DNAFile = [proj 'fasta/hm_h3n2_ny_ntcds_simple.fas'];
outFile = [proj 'beast/hm_h3n2_ny_dna_beast_1993_2006.fas'];
%generate_fasta_for_beast(DNAFile, merged_dat, outFile, 1993.601, 2006.6);
%disp 'BEAST input generated';