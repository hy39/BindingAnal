%Purpose: change the FASTA header
%Parameters: None
%Input files: FASTA files downloaded from Influenza Viruses Database
%Add the time information in FASTA header
%Output: 
%   FASTA file
function [] = generatet_fasta_for_beast_full(DNAFile, merged_dat, outFile, year_start, year_end)

    % metadata
    metadata = merged_dat;
    %nucleotide
    file_name = DNAFile; %'h1n1/hm_h1n1_noram_ntcds_simple1.fas';


    [Header, Sequence] = fastaread(['dat/' file_name]);
    New_Header = {};
    New_Sequence = {};
    
    %retrieve metadata
    load(['dat/' metadata]);
    if ~exist('gbacc_all')
        gbacc_all = cellstr(gbacc);
    end
    
    % I've changed the iso_date_num variable. The code only work for h3n2 ny Jul 10,2013 
    iso_date_number = iso_date_num.iso_date_number; %number in decimal format
    iso_date = iso_date_num.iso_date;               %original date format
    %ages = ages;
    ngs = Viruses(:,4);
    netcharges = Viruses(:,3);


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
    if (~isempty(id) & iso_date_number(id)>year_start & iso_date_number(id)<=year_end)

      %simple header
      %New_Header(length(New_Header)+1) = {[char(gbacc_all(id,2)) '_' num2str(isodate)]};
      %{protein accession, gene accession, isolation date, #NGS, netcharge} 
      
      %full header
      %New_Header(length(New_Header)+1) = {[char(Header(i)) '_' char(gbacc_all(id,2)) '_' num2str(iso_date(id)) '_' num2str(ngs(id)) '_' num2str(netcharges(id))]};
      genbank = '';
      %if isempty(gbacc_all) |  length(gbacc_all(1,:))==1
         New_Header(length(New_Header)+1) = {[char(Header(i)) '_' iso_date(id,:) '_' num2str(ngs(id)) '_' num2str(netcharges(id))]};
      %else    
      %   genbank = char(gbacc_all(id,2));
      %   New_Header(length(New_Header)+1) = {[char(Header(i)) '_' genbank '_' iso_date(id,:) '_' num2str(ngs(id)) '_' num2str(netcharges(id))]};
      %end
      New_Sequence(length(New_Sequence)+1) = Sequence(i);
    end
end

fastawrite(['dat/' outFile], New_Header, New_Sequence);
end


