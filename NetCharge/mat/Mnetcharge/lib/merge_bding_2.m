%Merge netcharge and binding affinity information.
%Written on October 23, 2012

function [] = merge_bding_2(bding_file, merged_bding_out)
Viruses = []; %(AGE, ISO_DATE, NETCHARGE, NGS)
header = [];
net_charge = [];
both_charge = [];

%Input parameters:
if ~exist('bding_file')
    %bding_file = 'h1n1/hm_h1n1_noram_charge.mat';
    disp 'error';
end

if ~exist('merged_bding_out')
    %merged_bding_out = 'h1n1/hm_h1n1_merged_data';
    disp 'error';  
end

%Retrieve binding data
fid = fopen(['dat/' bding_file]);
binding_info = textscan(fid, '%s %f %f %f', 'Delimiter', ',');


% Merge for every time stamped data
load(['dat/' merged_bding_out]);
if ~exist('gbacc_all')
 gbacc_all = cellstr(gbacc);
end

for i=1:length(gbacc(:,1))
 indx = find(strcmp(gbacc(i,:),cellstr(binding_info{1,1}))==1);
 %Viruses(i,1) = ages;
 %Viruses(i,2) = iso_date_num;
 %Viruses(i,3) = net_charge;
 %Viruses(i,4) = ngs;
 %Viruses(i,5) = both_charge;
 %Viruses(i,6) = infectionK;
 if isempty(indx)
     continue;
 end
 Viruses(i,7) = binding_info{1,2}(indx);
 Viruses(i,8) = binding_info{1,3}(indx);
 Viruses(i,9) = binding_info{1,4}(indx);
end


% Merge for data with aga information
save(['dat/' merged_bding_out], 'Viruses', 'gbacc', 'gbacc_all', 'ages', 'iso_date_num', 'binding_info');
end

