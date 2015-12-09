function [ ] = fasta2mat( filename )
%FASTA2MAT Summary of this function goes here
%   Detailed explanation goes here
% Transfer fasta sequences into mat variables
% >CY006211,80977927,A/Memphis/1/1968,1968//,4 (HA)

ID_date = {};
ID_year = [];
ID_num = [];
ID_nodate = [];
n_sequences = 0;
sequences = [];
header = {};

[Header Sequences] = fastaread(filename);
seq_length = length(char(Sequences(1)));
for i=1:length(Header)
  if length(char(Sequences(i))) ~= seq_length
    continue;
  end
  header_cell = regexp(char(Header(i)),',','split');
  access =  regexprep(char(header_cell(1)),'cds:','');
  gi = char(header_cell(2));
  strain = char(header_cell(3));;
  iso_date = char(header_cell(4));;
  d = regexp(iso_date,'/','split');
  yyyy = char(d(1));
  mm = char(d(2));
  dd = char(d(3));
  if isempty(mm)
      mm = '07';
      dd = '01';
      iso_date = [yyyy '/' mm '/' dd];
      ID_nodate(end+1,:) = str2num(gi);
  elseif isempty(dd)
      dd = '01';
      iso_date = [yyyy '/' mm '/' dd];
      ID_nodate(end+1,:) = str2num(gi);
  end
      
      
  segment = char(header_cell(5));;

  ID_num(end+1,:) = str2num(gi);
  ID_date(end+1,:) = {char(iso_date)};
  ID_year(end+1,:) = str2num(yyyy);
  sequences(end+1,:) = char(Sequences(i));
  header(end+1) = {[access '_' gi '_' iso_date]};
  n_sequences = n_sequences+1;
end

output = [filename(1:regexp(filename, '\.')-1) '.mat']; 
save(output, 'ID_date','ID_year','ID_num','n_sequences','sequences','header','ID_nodate');
end

