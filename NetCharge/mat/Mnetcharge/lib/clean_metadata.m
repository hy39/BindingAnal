%Merge netcharge and age/date with nonredundant information.
%Check Lab Note: September 7, 2012
%v2: Feb 4, 2013
%remove the redundant sets

function [gbacc ages iso_date] = clean_metadata(metadata, redun_set)


Filename1 = metadata;
% Retrieve metadata
fid = fopen(['dat/' Filename1]);
patient_info = textscan(fid, '%s %d %s', 'Delimiter', ',');
gbacc = char(patient_info{1,1});
ages = double(patient_info{1,2});
iso_date = char(patient_info{1,3});
cstr = cellstr(gbacc);
idx = find(strcmp(cstr, 'ADM07504'));

Filename2 = redun_set;
% Retrieve metadata
fid = fopen(['dat/' Filename2]);
patient_info = textscan(fid, '%s %s', 'Delimiter', ',');
gbacc_remove = char(patient_info{1,1});
strain_name = char(patient_info{1,2});

for i=1:length(gbacc_remove(:,1))
  idx = find(strcmp(cstr, gbacc_remove(i,:)));
  if (idx>0)
      gbacc(idx,:) = [];
      ages(idx,:) = [];
      iso_date(idx,:) = [];
  end
end
C = [cellstr(gbacc) cellstr(num2str(ages)) cellstr(iso_date)];

end

