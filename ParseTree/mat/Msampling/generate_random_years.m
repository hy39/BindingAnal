%Randomly select a subset of sequence from a larger set.
%filename: enter file name you want to open
%n_sequences: set size of current file
%siz: set desired size of new file

function generate_random_years(filename, n, s_year, e_year)
global Flu_Date;
global Flu_Number;
global s_year;
%s_year = 1980;
%e_year = 2008;
maxno = 15;
Flu_Date = [];

%% Read flu data and store into Flu_Date Object
%filename = 'B_1041.mat';
%filename = 'hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds.mat';
%filename2 = 'H1N1_1032_aligned_subset.mat';
%filename3 = 'H3N2_987_aligned_subset.mat';
%filenames = strvcat(filename1,filename2,filename3);
%for i1 = 1:3
load(filename,'-mat');
Flu_Date.ID_year = ID_year;
%end

Flu_Number = SetMaxSize(Flu_Date, s_year, e_year, maxno);
%FindMaxSize(2000);

%meta data for random years files
file_name = 'id';
load(['dat/years/' file_name '.mat'],'-mat');
fid = fid + 1;
data = [];

%input total number of sequences
n_desired = n;

%randomLocs = randperm(n_sequences);
y = s_year + (e_year-s_year).*rand(n_desired,1);
y = round(y);
%if (find(y==1993) > 0)
%  y(find(y==1993)) = round(1977 + (2007-s_year).*rand(1,1));
%end
sy = sort(y);

%check whether the sample size over the limit for each year 
for j=1:length(sy) % need to modify this line. sy not necessarily to be continue 
    while(length(find(y==sy(j)))>FindMaxSize(sy(j), Flu_Number))
        locs = find(y==sy(j));
        disp('test');
        %if the sample size over the limit, move the sampled year to the next year
        y(locs(end)) = sy(j)+1;%increase 1 for the selected year
        sy=sort(y);%sort the year again
    end
end





%save randseq
output = ['dat/years/ys_' int2str(n) '_' int2str(fid) '.mat']; 
save(output,'Flu_Date','Flu_Number','s_year','e_year','fid','filename','n_desired','y','sy');
yr_data(fid).id = fid;
yr_data(fid).fname = output;
save('dat/years/id', 'fid', 'output','yr_data');
clear all;
end


function  [m] = FindMaxSize(yr, Flu_Number)
m = Flu_Number(find(Flu_Number(:,1) == yr),end);
end



%%To be done
%%After random generate years, calculate the maximize strains size that
%%included in all of 3 types of virues in the selected year
function [Flu_Number] = SetMaxSize(Flu_Date, s_year, e_year, maxn)
%Flu_Number = zeros(1,4);
cnt = 0;
for ii2 = s_year:e_year %2008
  cnt = cnt + 1;
  Flu_Number(cnt,1) = ii2;
  for fileid = 1:length(Flu_Date)
    Flu_Number(cnt,fileid+1) = length(find(Flu_Date(fileid).ID_year == ii2));
  end
  Flu_Number(cnt,length(Flu_Date)+2) = min([Flu_Number(cnt,2:length(Flu_Date)+1) maxn]);
%  Flu_Number(cnt,length(Flu_Date)+2) = 1;
end
end

