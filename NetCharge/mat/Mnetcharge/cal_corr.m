% plot the netcharge vs ages
ranges = [1993.601 1994.6; 
          1994.601 1995.6;
          1995.601 1996.6; 
          1996.601 1997.6;
          1997.601 1998.6; 
          1998.601 1999.6;
          1999.601 2000.6; 
          2000.601 2001.6;
          2001.601 2002.6;
          2002.601 2003.6; 
          2003.601 2004.6;
          2004.601 2005.6; 
          ];
%ranges2 = [1993.601 2006.6];

%%%ranges2 = [2003.601 2004.6]; %consistent 
%%%ranges2 = [2001.601 2002.6]; %consistent     
cor = [];
for i=1:length(ranges(:,1))
    ranges2 = ranges(i,:);
    cor(i,1:2) = ranges2; 
    %ranges2 = [2004.601 2005.6];   
    load('merged_data.mat');



    ngs_index = 4;
    charge_index = 3;
    index = 1; %1)age, 2)isolation date
    year_s = 1;
    % #NGS=8
    TF = find(Viruses(:,2)>ranges2(1,1) & Viruses(:,2)<ranges2(1,2) &  Viruses(:,1)~=0 & Viruses(:,1)<70);
    viruses = Viruses(TF,:);

    %hist(viruses(:,1),10);figure(gcf);
    ranges2
    r = corrcoef(viruses(:,1), viruses(:,3)); %correlation age vs netcharge
    if length(r)==2
        cor(i,3) = r(1,2);
    else 
        cor(i,3) = 1;
    end
    r = corrcoef(viruses(:,2), viruses(:,3)); %correlation isolated date vs netcharge
    if length(r)==2
        cor(i,4) = r(1,2);
    else
        cor(i,4) = 1;    
    end
    r = corrcoef(viruses(:,1), viruses(:,4)); %correlation age vs #NGS
    if length(r)==2
        cor(i,5) = r(1,2);
    else
        cor(i,5) = 1;    
    end
    
end
save('cirrelation', 'cor');