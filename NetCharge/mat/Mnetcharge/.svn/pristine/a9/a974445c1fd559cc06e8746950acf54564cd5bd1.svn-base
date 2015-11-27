%Merge binding avidity score into mat file
%Feb 22, 2012

function [Viruses gbacc] = merge_binding(binding_dat, merged_dat, merged_out)
Viruses = []; %(AGE, ISO_DATE, NETCHARGE, NGS)
gbacc_all;

binding_gb = []; 
binding_scr = [];
structure_scr = [];
sum_scr = [];



% Retrieve netcharge
dat = load(['dat/' binding_dat]);
header = dat.header_txt;
net_charge = dat.charge_txt(:,1);
both_charge = dat.charge_txt(:,4);
ngs = dat.charge_txt(:,5);

% Retrieve metadata
fid = fopen(['dat/' Filename2]);
patient_info = textscan(fid, '%s %d %s', 'Delimiter', ',');
gbacc = char(patient_info{1,1});
ages = double(patient_info{1,2});
iso_date = char(patient_info{1,3});
iso_date_num = [];
%iso_date_str = [];

% convert date_string into numeric date
for i=1:length(iso_date(:,1))
date_str = strread(iso_date(i,:), '%d', 'delimiter', '/');
if (length(date_str) == 3)
    yy = date_str(1,1);
    mm = date_str(2,1);
    dd = date_str(3,1);
  else 
    yy = date_str(1,1);
    mm = 0;
    dd = 0;
end
date_num = yy+(datenum(yy,mm,dd)-datenum(yy,1,1))./(datenum(yy,12,31)-datenum(yy,1,1));
iso_date_num(i,1) = round(10^3*date_num)/10^3;
%iso_date_str(i,:) = char(iso_date_num(i,1));
end
%n2 = datenum(,'yyyy/mm/dd')

% Merge for every time stamped data
for i=1:length(gbacc(:,1))
TF = strcmp(gbacc(i,:),cellstr(header));
Viruses(i,1) = ages(i);
Viruses(i,2) = iso_date_num(i);
Viruses(i,3) = net_charge(TF);
Viruses(i,4) = ngs(TF);
Viruses(i,5) = both_charge(TF);
end

% Merge for data with aga information


save(['dat/' outFile], 'Viruses', 'gbacc', 'ages', 'iso_date_num');
end

