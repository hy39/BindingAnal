%Merge netcharge and age/date information.
%Check Lab Note: September 7, 2012

function [Viruses gbacc] = merge_metadata(charge_dat, metadata, merged_out)
Viruses = []; %(AGE, ISO_DATE, NETCHARGE, NGS)
header = [];
net_charge = [];
both_charge = [];

%Input parameters:
if ~exist('charge_dat')
    disp 'error: not enough parameters for charged data.';
else
    Filename1 = charge_dat;
end

if ~exist('metadata')
    disp 'error: not enough parameters for metadata.';
else
    Filename2 = metadata;
end

if ~exist('merged_out')
    outFile = 'h1n1/hm_h1n1_merged_data';
else
    outFile = merged_out;
end

% Retrieve netcharge
dat = load(['dat/' Filename1]);
header = dat.header_txt;
net_charge = dat.charge_txt(:,1);
both_charge = dat.charge_txt(:,4);
ngs = dat.charge_txt(:,5);

% Retrieve metadata
fid = fopen(['dat/' Filename2]);
patient_info = textscan(fid, '%s %f %s', 'Delimiter', ',');
gbacc = char(patient_info{1,1});
ages = double(patient_info{1,2});
iso_date = char(patient_info{1,3});
iso_date_number = [];
iso_date_number_int = [];
%iso_date_str = [];

% convert date_string into numeric date
for i=1:length(iso_date(:,1))
    date_str = strread(iso_date(i,:), '%d', 'delimiter', '/');
    iso_date_str = '';
    if (length(date_str) == 3)
        yy = date_str(1,1);
        mm = date_str(2,1);
        dd = date_str(3,1);
        iso_date_str = iso_date(i,:);
    else % if there are no month and day information, assume Jul 1.
        yy = date_str(1,1);
        mm = 7;
        dd = 1;
        iso_date_str = [num2str(yy) '/07/01'];
    end
    date_num = yy+(datenum(yy,mm,dd)-datenum(yy,1,1))./(datenum(yy,12,31)-datenum(yy,1,1));
    iso_date_number(i,1) = round(10^3*date_num)/10^3;
    %iso_date_number_int(i,1) = datenum(iso_date(i,:),'yyyy/mm/dd');
    iso_date_number_int(i,1) = datenum(iso_date_str,'yyyy/mm/dd');
end

% Merge for every time stamped data
for i=1:length(gbacc(:,1))
TF = strcmp(gbacc(i,:),cellstr(header));
if isempty(find(TF==1))
  continue;
end
Viruses(i,1) = ages(i);
Viruses(i,2) = iso_date_number(i);
Viruses(i,3) = net_charge(TF);
Viruses(i,4) = ngs(TF);
Viruses(i,5) = both_charge(TF);
end

% Merge for data with aga information
iso_date_num = struct;
iso_date_num.iso_date = iso_date;
iso_date_num.iso_date_number = iso_date_number;
iso_date_num.iso_date_number_int = iso_date_number_int;

save(['dat/' outFile], 'Viruses', 'gbacc', 'ages', 'iso_date_num');
end

