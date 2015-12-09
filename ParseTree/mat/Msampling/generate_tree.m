%Randomly select a subset of sequence from a larger set.
%filename: enter file name you want to open
%n_sequences: set size of current file
%siz: set desired size of new file

function generate_tree(filename)


load(filename,'-mat');
%clear data;
clear siz;
data = [];
storedData = [];

%n_desired = num;
n_desired = length(ID_num);
for i1 = 1:n_desired
            data(i1).Header = char(header(i1));
            data(i1).GI = int2str(ID_num(i1)); % need to use dataSiz(1,2) to obtain the size
            data(i1).Year = ID_year(i1);
            data(i1).Sequence = sequences(i1,:);
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

% create output directory
out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7)]
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end

% write fasta file
spos = regexp(filename, '\/');
file_type = filename(spos(end):regexp(filename, '\.')-1);
file_sample = ['whole_' num2str(n_desired)];
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

