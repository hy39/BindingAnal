fileID = fopen('taxa.txt','w');
fprintf(fileID,'%s', taxa(1,189).annotation)
fclose(fileID);


fileID = fopen('h3n2_ny_binding_n1371.csv');
A = fscanf(fileID,'%d %f %f %f');
fclose(fileID);


M = csvread('h3n2_ny_binding_n1371.csv');
binding_txt = M(:,1:2);
save('hm_h3n2_ny_binding_1993.mat','binding_txt');