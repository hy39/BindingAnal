%Purpose: change the FASTA header
%Parameters: None
%Input files: FASTA files downloaded from Influenza Viruses Database
%Add the time information in FASTA header
%Output: 
%   FASTA file
function [] = h1n1_format_fasta_header(DNAFile, merged_dat, outFile, year_start, year_end)
%Human H3N2 New York Strains

if ~exist('merged_dat')
    metadata = 'hm_h1n1_merged_data.mat';
else
    metadata = merged_dat;
end
%protein
%file_name = 'hm_h3n2_flu_ny_any_simple.fas';
%out_file_name = 'hm_h3n2_ny_beast.fas';
%out_file_name = 'hm_h3n2_ny_age_beast_1993_2005.fas'; % ages>0

%nucleotide
if ~exist('DNAFile')
    file_name = 'h1n1/hm_h1n1_noram_ntcds_simple1.fas';
else
    file_name = DNAFile;
end

if ~exist('outFile')
    out_file_name = 'h1n1/hm_h1n1_noram_age_dna_beast_2000_2008.fas'; % ages>0
else
    out_file_name = outFile;
end

[Header, Sequence] = fastaread(['dat/' file_name]);
New_Header = {};
New_Sequence = {};

%retrieve metadata
dat = load(['dat/' metadata]);
gbacc_all = dat.gbacc_all;
iso_date = dat.iso_date_num;
ages = dat.ages;

siz = length(Header);
for i=1:siz
    %%Formatting
    id = find(strcmp(gbacc_all(:,1),Header(i))==1);
    
    if isempty(id)
      continue;
    end
    
    %length(New_Header)+1
    isodate = [];
    %if (~isempty(id) & ages(id)~=0 & iso_date(id)>2000 & iso_date(id)<=2008)
    if (~isempty(id) & iso_date(id)>year_start & iso_date(id)<=year_end)
      isodate = iso_date(id);
      %New_Header(length(New_Header)+1) = {[char(Header(i)) '_' num2str(isodate)]};
      New_Header(length(New_Header)+1) = {[char(Header(i)) '_' char(gbacc_all(id,2)) '_' num2str(isodate)]};
      New_Sequence(length(New_Sequence)+1) = Sequence(i);
    end
end

fastawrite(['dat/' out_file_name], New_Header, New_Sequence);
end


