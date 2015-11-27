%Purpose: calculate the number of each amino acid occured in given sequences
%Parameters: None
%Input files: FASTA files downloaded from Influenza Viruses Database
%Output: 
%   1) Text file: {header, net charge, positive charge, negative chage, both} 
%   1) aa_array: cell array of #amino acids for each fasta sequence
%   2) aa_charge: isolation year and net charges for each fasta sequence
%   3) aa_charge_sum: year, net charges 
function [header_txt charge_txt] = main_aa_dist(infile, outfile, subtype)

    file_name = infile;
    out_file_name = outfile;
    
%--------------------------------------------------------------------------
[Header, Sequence] = fastaread(['dat/' file_name]);
siz = length(Header);
aa_array = cell(siz,1);
%s = struct([])
aa_charge = [];
aa_charge_sum = [];
for i=1:siz
    aa_seq = char(Sequence(i))
    aa = get_aa_num(aa_seq);
    
    %specify the region in protein sequences
    if strcmp(subtype,'H3N2')
        aa = get_aa_num(aa_seq(17:345)); %HAR1. For H3N2
    elseif strcmp(subtype,'H1N1')
        aa = get_aa_num(aa_seq(1:323)); %H1N1 clean. For H1N1
    end
    %aa = get_aa_num(aa_seq(1:318)); %H3N2 clean
    
    aa_array{i} = aa;
    aa_charge(i,1) = aa.R + aa.H + aa.K - aa.D - aa.E; %Net charge
    aa_charge(i,2) = aa.R + aa.H + aa.K; % positive charge
    aa_charge(i,3) = aa.D + aa.E;        % negative charge
    aa_charge(i,4) = aa_charge(i,2) + aa_charge(i,3); % #both
    aa_charge(i,5) = aa.NGS; % #NGS
    %--
    ngs(i).header = char(Header(i));
    ngs(i).number = aa.NGS;
    ngs(i).location = aa.locs;
end

    header_txt = char(Header');
    charge_txt(:,1) = aa_charge(:,1);
    charge_txt(:,2) = aa_charge(:,2);
    charge_txt(:,3) = aa_charge(:,3);
    charge_txt(:,4) = aa_charge(:,4);
    charge_txt(:,5) = aa_charge(:,5);
    
    save(['dat/' out_file_name], 'header_txt', 'charge_txt', 'ngs');
    
    %save to text file
    fileID = fopen(['dat/' out_file_name '.csv'],'w');
    for i=1:length(header_txt)
    fprintf(fileID,'%8s,%2d,%2d,%2d,%2d\n',header_txt(i,:), charge_txt(i,1), charge_txt(i,2), charge_txt(i,3), charge_txt(i,4));
    end
    fclose(fileID);
end


