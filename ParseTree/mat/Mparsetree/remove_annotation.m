function [tr] = remove_annotation( )
file = 'dat/ancestor_20130702/hm_h3n2_ny_dna_beast_1993_2006_test.traits.nw.mcc.(time).trees';
ff = fopen(file,'r');
[tr, k] = fscanf(ff,'%c');

% remove annotation
while ~ isempty(regexp(tr, '['))
    pos = regexp(tr, '[');
    pos_s = pos(1);
    pos = regexp(tr, ']');
    pos_e = pos(1);
    tr(pos_s:pos_e) = [];
end

% remove branch lengths
tr = regexprep(tr,':-?\d+.\d+,',',');
tr = regexprep(tr,':-?\d+.\d+)',')');

fileID = fopen('tree.txt','w');
fprintf(fileID,'%s',tr);
fclose(fileID);

