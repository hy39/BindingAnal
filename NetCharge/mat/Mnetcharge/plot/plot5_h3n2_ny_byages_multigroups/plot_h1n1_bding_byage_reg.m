function [] = plot_h1n1_bding_byage_reg( age_min, age_max, saveoutput )
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
ranges2 = [1993.601 2006.6]; %H1N1: 1993-2006 
%ranges2 = [2005.601 2008.6];   
load('../../dat/h1n1_noram_1995/hm_h1n1_noram_merged_data_noinsert.mat');
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
    
    viruses_ngs = viruses(find(viruses(:,ngs_index)>=7 & viruses(:,ngs_index)<=11 & viruses(:,ageid)<=ages_lim & viruses(:,ageid)>=ages_min),:);
    x_tot = viruses_ngs(:,index);
    y_tot = viruses_ngs(:,Y_IDX);
    x = x_tot;
    
    hold on;
    color_code = strvcat('G', 'B', 'R', 'K');
    for ngs_no=7:7
        ngs_no
        dm1 = find(viruses_ngs(:,ngs_index)==ngs_no);
        loc_x1 = find(viruses_ngs(:,ngs_index)==ngs_no);
        x1 = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no),index);
        y1 = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no),Y_IDX);
        %plot(x1, y1, [color_code(ngs_no-6) 'o'], 'MarkerSize', 2, 'MarkerFaceColor',[color_code(ngs_no-6)]);
   
        % detect outliers
        reg_stats = regstats(y1,x1);
        potential_outlier = reg_stats.cookd > 4/length(x1);
        %scatter(x1(potential_outlier),y1(potential_outlier), 'rx', 'LineWidth', 0.5)
        % drop outliers
        x1(potential_outlier) = [];
        y1(potential_outlier)= [];
        viruses_ngs(loc_x1(potential_outlier),10) = 1;
        
        plot(x1, y1, [color_code(ngs_no-6) 'o'], 'MarkerSize', 3, 'MarkerFaceColor',[color_code(ngs_no-6)]);
        X = [ones(size(x1)) x1]; 
        [b,bint,r,rint,stats] = regress(y1,X);
        x1fit = Xmin:ages_lim;
        YFIT = b(1) + b(2)*x1fit;
        plot(x1fit,YFIT, [color_code(ngs_no-6) '--'], 'LineWidth', 2);
        [slope pv se T] = getPfromTtest(X,y1);
        slope = round(slope.*1000.*100)./100;
        se = round(se.*1000.*100)./100;
        j1 = length(pvalue)+1;
        %pvalue(j1).metadata = 'season start, season end, slope, upper bound, lower bound, R2, F, p value';
        pvalue(j1).metadata = 'season start, season end, slope, T score, p value, standard error, R2';
        pvalue(j1).season = ranges2;
        pvalue(j1).ngs = ngs_no;
        pvalue(j1).slope = b(2);
        pvalue(j1).bint = bint;
        pvalue(j1).stats = stats;
        pvalue(j1).p = stats(3);
        pvalue(j1).cookd = reg_stats.cookd;
        pvalue(j1).n = length(X);
        %ptable(j1,:) = [pvalue(j1).season pvalue(j1).ngs pvalue(j1).slope bint(2,1) bint(2,2) pvalue(j1).stats(1) pvalue(j1).stats(2) pvalue(j1).stats(3)];
        ptable(j1,:) = [pvalue(j1).season pvalue(j1).ngs slope T pv se stats(1)];
    end
    if exist('saveoutput')
        if strcmp(saveoutput,'s')
            save(['h1n1_ny_bding_byage' num2str(ages_min) '-' num2str(ages_lim) '_reg_pvalue.mat'], 'pvalue', 'ptable', 'ages_min', 'ages_lim');
        end
    end
    pvalue = [];
    ptable = [];
    
    %create dummy variable
    viruses_ngs_noutlier = viruses_ngs(find(viruses_ngs(:,10)~=1),:);
    ngs = nominal(viruses_ngs_noutlier(:,4));
    dv = dummyvar(ngs);
    dv(:,1) = 1;
    x_tot = viruses_ngs_noutlier(:,index);
    y_tot = viruses_ngs_noutlier(:,Y_IDX);
    
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
    %Y1FIT = b(1) + b(end)*x1fit;
    %plot(x1fit,Y1FIT, 'G-','LineWidth',1.5);
    %Y1FIT = b(1)+b(2) + b(end)*x1fit;
    %plot(x1fit,Y1FIT, 'B-','LineWidth',1.5);
    %Y1FIT = b(1)+b(3) + b(end)*x1fit;
    %plot(x1fit,Y1FIT, 'R-','LineWidth',1.5);
    %Y1FIT = b(1)+b(4) + b(end)*x1fit;
    %plot(x1fit,Y1FIT, 'K-','LineWidth',1.5);
    
    
    
    
    xlabel('host age(in yrs)');
    ylabel('binding score');
    xlim([0 round(Xlim./10)*10]);
    
    %if exist('saveoutput')
    %    if strcmp(saveoutput,'s')
    %        save(['h1n1_ny_bding_byage' num2str(ages_min) '-' num2str(ages_lim) '_dummy_pvalue.mat'], 'pvalue', 'ptable', 'ages_min', 'ages_lim');
    %    end
    %end
end  

