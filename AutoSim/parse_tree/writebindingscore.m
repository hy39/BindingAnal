proj = 'ancestor_20130713';
M = csvread(['dat/' proj '/h3n2_ny_binding_n1371.csv']);
binding_txt = M(:,[1 4]);
save(['dat/' proj '/hm_h3n2_ny_binding_1993.mat'],'binding_txt');