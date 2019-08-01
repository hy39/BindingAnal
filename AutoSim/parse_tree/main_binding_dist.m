function [binding_txt] = main_binding_dist(inHeader, inBinding, outfile)
% Example: main_binding_dist('ancestor_20130713/hm_h3n2_ny_charge_1993.mat','ancestor_20130713/h3n2_ny_binding_n1371.csv','ancestor_20130713/hm_h3n2_ny_binding_1993.mat')

    load(['dat/' inHeader]);
    fid = fopen(['dat/' inBinding],'r');
    bd = textscan(fid,'%s %f %f %f','delimiter', ',')
    binding_txt = [];
    for i=1:length(header_txt(:,1))
      id = str2num(header_txt(i,:));
      binding_txt(id,:) = [id bd{4}(strcmp(bd{1},num2str(id)))];
    end
    save(['dat/' outfile], 'binding_txt');
    
    %save to text file
    %fileID = fopen(['dat/' out_file_name '.csv'],'w');
    %for i=1:length(header_txt)
    %fprintf(fileID,'%8s,%2d,%2d,%2d,%2d\n',header_txt(i,:), charge_txt(i,1), charge_txt(i,2), charge_txt(i,3), charge_txt(i,4));
    %end
    %fclose(fileID);
end


