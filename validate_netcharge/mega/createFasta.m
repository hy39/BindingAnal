
T1 = readtable('test2.csv');
H1seq = '----DTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCRLKGIAPLQLGKCNIAGWLLGNPECDPLLPVRSWSYIVETPNSENGICYPGDFIDYEELREQLSSVSSFERFEIFPKESSWPNHNTNGVTAACSHEGKSSFYRNLLWLTEKEGSYPKLKNSYVNKKGKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVRDQAGRMNYYWTLLKPGDTIIFEANGNLIAPMYAFALRRGFGSGIITSNASMHECNTKCQTPLGAINSSLPYQNIHPVTIGECPKYVRSAKLRMVTGLRNIPAR----';

fasta = struct;
fasta(1).Header = ['>H1_1934', char(10)];
fasta(1).Sequence = str_new;
for i=1:length(T1.POS)
  h = T1{i,'AAPOS'};
  pos = T1{i,'POS'};
  wild = T1{i,'WILD'};
  mut = T1{i,'MUTANT'};
  str_new = '';
  if strcmp(H1seq(pos),wild)
        str_new = H1seq;
        str_new(pos) = char(mut);
        if ~strcmp(T1{i,'POS2'},'NA')
	    	pos2 = T1{i,'POS2'};
            pos2 = str2num(pos2{1});
            if isnumeric(pos2)
            wild2 = T1{i,'WILD2'};
            mut2 = T1{i,'MUTANT2'};
            str_new(pos2) = char(mut2);
            end
        end
        fasta(i+1).Header = ['>', char(h), char(10)];
        fasta(i+1).Sequence = str_new;
  end    
end
fastawrite('h1_netcharge.fas',fasta);