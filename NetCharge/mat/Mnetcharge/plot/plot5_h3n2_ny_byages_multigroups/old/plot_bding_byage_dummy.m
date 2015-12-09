function [] = plot_bding_byage_dummy (age_min, age_max, saveoutput)
% plot the binding vs netcharge
% Feb 28, 2013

addpath(genpath('..'));
addpath(genpath('../../lib'));
ranges = [
          1993.601 1994.6; 
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
          ]
%ranges2 = [1993.601 2005.6];
dominance = [9 0;
             9 0;
             8 0;
             8 0;
             10 0;
             11 10;
             10 0;
             11 0;
             10 0;
             11 10;
             10 11;
             11 0];
%%%ranges2 = [2003.601 2004.6]; %consistent 
%%%ranges2 = [2001.601 2002.6]; %consistent     
ranges2 = [1993.601 2005.6];   
load('../../dat/h3n2_ny/hm_h3n2_ny_merged_data.mat');
i = 1;


y_lim = [30 45];

pvalue = [];
ptable = [];


    ngs_index = 4;
    charge_index = 3; %netcharge
    %charge_index = 5; %total charged amino acids
    ageid = 1;
    kid = 6;
    %index = charge_index; % 1:age, 2:isolation date, 6:times of previouis infection, 3:netcharge 
    index = ageid;
    bding_charge_index = 7;
    bding_total_index = 9;
    Y_IDX = bding_total_index;
    %Y_IDX = charge_index;
    min_no = 12;
    
    % #NGS=8
    TF = find(Viruses(:,2)>ranges2(1,1) & Viruses(:,2)<ranges2(1,2) &  Viruses(:,1)~=0);
    viruses = Viruses(TF,:);

    
    dm_ngs = dominance(i,:);
    charge_lim = 20;
    ages_min = age_min;
    ages_lim = age_max;
    min_no = 12;
    Xlim = ages_lim;
    Xmin = 1;
    
   
    hold on;
    
    viruses_ngs = viruses(find(viruses(:,ngs_index)>=8 & viruses(:,ngs_index)<=11 & viruses(:,ageid)<=ages_lim & viruses(:,ageid)>=ages_min),:);
    x_tot = viruses_ngs(:,index);
    y_tot = viruses_ngs(:,Y_IDX);
    x = x_tot;
    
     color_code = strvcat('G', 'B', 'R', 'K');
    for ngs_no=8:11
    ngs_no
    dm1 = find(viruses_ngs(:,ngs_index)==ngs_no);
    dmI1 = zeros(size(x));
    dmI1(dm1)=1;
    x1 = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no),index);
    y1 = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no),Y_IDX);
    plot(x1, y1, [color_code(ngs_no-7) 'o'], 'MarkerSize', 2, 'MarkerFaceColor',[color_code(ngs_no-7)]);
    end
    

    
    
    %create dummy variable
    ngs = nominal(viruses_ngs(:,4));
    dv = dummyvar(ngs);
    dv(:,1) = 1;
    X = [dv x_tot];
    [b,bint,r,rint,stats] = regress(y_tot,X);
    %b(1):b(4) interception b(5) slope
    [slope pv se T] = getPfromTtestdv(X,y_tot,dv);
       
    
    j1 = length(pvalue)+1;
    pvalue(j1).metadata = 'season start, season end, slope, T score, p value, standard error, R2';
    pvalue(j1).season = ranges2;
    pvalue(j1).ngs = ngs_no;
    pvalue(j1).slope = slope;
    ptable(j1,:) = [pvalue(j1).season ngs_no slope T pv se stats(1)];
    
    
    x1fit = Xmin:Xlim;
    
    Y1FIT = b(1) + b(end)*x1fit;
    plot(x1fit,Y1FIT, 'R-');
    Y1FIT = b(1)+b(2) + b(end)*x1fit;
    plot(x1fit,Y1FIT, 'G-');
    Y1FIT = b(1)+b(3) + b(end)*x1fit;
    plot(x1fit,Y1FIT, 'B-');
    Y1FIT = b(1)+b(4) + b(end)*x1fit;
    plot(x1fit,Y1FIT, 'K-');
    
    xlabel('ages');
    ylabel('binding score');
    xlim([Xmin Xlim]);
    if exist('saveoutput')
        if strcmp(saveoutput,'s')
            save(['h3n2_ny_bding_byage' num2str(ages_min) '-' num2str(ages_lim) '_dummy_pvalue.mat'], 'pvalue', 'ptable', 'age_min', 'age_max');
        end
    end    
end
    
  

