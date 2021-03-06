%main script to analyze netcharge in viral phylogenies
%1)calculate netcharge -> charge.mat
%2)merge metadata with netcharge -> merged_data.mat
%3)add gene accession to each protein accession
%4)create DNA FASTA files for BEAST input
%Oct 2, 2012
%New Zealand 1993.6-2006.6
p = path
path(p,'lib/');

proj = 'h3n2_nz/';

%%%1)calculate charge distribution
%prot_fasta = [proj 'fasta/hm_h3n2_nz_any_simple.fas'];
%charge_out = [proj 'hm_h3n2_nz_charge'];
%%[header charge] = main_aa_dist(prot_fasta, charge_out, 'H3N2');
disp 'Scanning amino acids';

%%%2)merge metadata with netcharge -> merged_data.mat
charge_dat = [proj 'hm_h3n2_nz_charge'];
metadata = [proj 'hm_h3n2_nz_any_ay.csv'];
merged_out = [proj 'hm_h3n2_nz_merged_data'];
%[v gbacc] = merge_metadata(charge_dat, metadata, merged_out);
disp 'Charged data generated';

%%%3)add gene accession to each protein accession
merged_dat = [proj 'hm_h3n2_nz_merged_data'];
accession_table_out = [proj 'hm_h3n2_nz_any'];
%create_geneprot_table(accession_table_out, gbacc, merged_dat); %No need
%%convert_age2k(merged_dat);                                    %No need
disp 'Add GenBank gene accession';

%%%4)merge with predicted binding avidity data
merged_dat = [proj 'hm_h3n2_nz_merged_data'];
bding = [proj 'scr/h3n2_nz_binding.csv'];
%merge_bding_2(bding, merged_dat);
disp 'Binding data generated';

%%%5)create DNA FASTA files for BEAST input
DNAFile = [proj 'fasta/hm_h3n2_nz_ntcds_simple.fas'];
outFile = [proj 'beast/hm_h3n2_nz_dna_beast_1993_2006.fas'];
generate_fasta_for_beast(DNAFile, merged_dat, outFile, 1993.601, 2006.6);
disp 'BEAST input generated';