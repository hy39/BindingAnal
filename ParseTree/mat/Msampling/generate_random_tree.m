%Randomly select a subset of sequence from a larger set.
%filename: enter file name you want to open
%n_sequences: set size of current file
%siz: set desired size of new file

function generate_random_tree(filename, ry)
if (nargin == 1)
  %num = 1;
   num = '1';
end;
%filename = 'B_1041.mat';
%filename = 'B_1041_aligned_subset';
%filename = 'H1N1_1032_aligned_subset';
%filename = 'H3N2_987_aligned_subset';

load(filename,'-mat');
%clear data;
clear siz;
data = [];
storedData = [];


%open random generated years from dat/years.mat
filename_yrs = ry;
load(['dat/years/' filename_yrs],'-mat');

%if no random generated years, generate it now
%input total number of sequences
%randomLocs = randperm(n_sequences);

%n_desired = 25; %change to 300
%y = 1980 + (2008-1980).*rand(n_desired,1);
%y = round(y); %random sampled years


%find the locations for each random generated year
%need to pay more attention on some year with only few records.
%it will cause final data size less than desired number

for i1 = 1:n_desired
    locs = find(ID_year == y(i1));
    %random location of each strain's element
    randomLocs = locs(randperm(length(locs)));
    locs_siz = size(locs);
    j1 = 1;
    
    %%%%To be done: Check whether all of the records in the selected years
    %%%%            have been chosen 
    %%%%
    %while (isempty(randomLocs) || locs_siz(1) <5)
    while (isempty(randomLocs))
         %find closest year in ID_year
         locs = find(ID_year == y(i1)+j1);
         randomLocs = locs(randperm(length(locs)));
         locs_siz = size(locs);
         j1=j1+1;
    end
    
    %check if there is same element stored in Data already
    %skip this element if this is chosen before
    currentElement = ID_num(randomLocs(1));
    for i2 = 1:length(randomLocs)
        currentElement = ID_num(randomLocs(i2)); % choose a GI
        if (isempty(find(storedData == currentElement)))
            dataSiz = length(data);
            %data(dataSiz+1).Header = int2str(ID_num(randomLocs(i2)));
            data(dataSiz+1).Header = char(header(randomLocs(i2)));
            data(dataSiz+1).GI = int2str(ID_num(randomLocs(i2))); % need to use dataSiz(1,2) to obtain the size
            data(dataSiz+1).Year = ID_year(randomLocs(i2));
            data(dataSiz+1).Sequence = sequences(randomLocs(i2),:);
            
            storeSiz = length(storedData);
            %storedData(storeSiz(1,2)+1) = currentElement;
            storedData = [storedData currentElement];
            break;
        end
    end
end

%sort by year
for i = 1:n_desired-1
    for j2 = 1:n_desired-1
        if data(j2).Year > data(j2+1).Year
            header_tmp = data(j2+1).Header;
            gi_tmp = data(j2+1).GI;
            sequence_tmp = data(j2+1).Sequence;
            year_tmp = data(j2+1).Year;
            data(j2+1) = data(j2);
            data(j2).Header = header_tmp;
            data(j2).GI = gi_tmp;
            data(j2).Sequence = sequence_tmp;
            data(j2).Year = year_tmp;
        end
    end
end

%%change header to be sequential ID chronologically
%%if desired number is < 100, then there should always be 2 digits
%%%data = update_header(data, n_desired); %change to sequential order




% create output directory
out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7)]
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end

% write fasta file
spos = regexp(filename, '\/');
file_type = filename(spos(end):regexp(filename, '\.')-1);
file_sample = filename_yrs(1:regexp(filename_yrs, '\.')-1);
output_fasta = [out_dir '/' file_type '_' file_sample '.fas'];
fastawrite(output_fasta,data);
clear;
end

% header will be using serial number starting from 1
function [data] = update_header(data, n_desired)
if(n_desired >= 10 && n_desired < 100)
for i = 1:n_desired
            if(i<10)
                data(i).Header = ['0' int2str(i)];
            else
                data(i).Header = int2str(i);
            end
end
end
%if desired number is < 100, then there should always be 3 digits
if(n_desired >= 100 && n_desired <= 1000)
for i = 1:n_desired
            if(i<100)
                if(i<10)
                data(i).Header = ['00' int2str(i)];
                else
                data(i).Header = ['0' int2str(i)];
                end
            else
                data(i).Header = int2str(i);
            end
end
end

end


function check_stored()
         locs = find(ID_year == y(i1)+j);
         randomLocs = locs(randperm(size(locs)));
         locs_siz = size(locs);
end


% Try do write a function can list ID mapping for GI and Serial#

%function [data] = update_header(data, n_desired)
%if(n_desired >= 10 && n_desired < 100)
%for i = 1:n_desired
%            if(i<10)
%                data(i).Header = ['0' int2str(i)];
%            else
%                data(i).Header = int2str(i);
%            end
%end
