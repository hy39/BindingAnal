function [] = plot_h1n1_net_byage( age_min, age_max, saveoutput )

p = path
path(p,'..');
p = path
path(p,'../lib');
p = path
path(p,'../../lib');

% plot the netcharge vs ages
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

%%%ranges2 = [2003.601 2004.6]; %consistent 
%%%ranges2 = [2001.601 2002.6]; %consistent     
ranges2 = [1993.601 2006.6]; %H1N1: 1993-2006 
%ranges2 = [1995.001 2006.6]; %H1N1: 1995-2008   
%ranges2 = [2005.001 2008.6]; %H1N1: 1995-2008   
%load('../dat/h3n2_ny/hm_h3n2_ny_merged_data.mat');
load('../../dat/h1n1_noram_1995/hm_h1n1_noram_merged_data_ex5.mat'); %exclude 5 strains
i = 1;
pvalue = [];
ptable = [];

for i=1:length(ranges2(:,1))
    ranges2 = ranges2(i,:);
    ngs_index = 4;
    charge_index = 3; %netcharge
    %charge_index = 5; %total charged amino acids
    K_index = 6;
    index = 1; %1)age, 2)isolation date
    ageid = 1;
    % #NGS=8
    TF = find(Viruses(:,2)>ranges2(1,1) & Viruses(:,2)<ranges2(1,2) &  Viruses(:,1)~=0);
    viruses = Viruses(TF,:);
    Y_IDX = charge_index;
    
    ages_min = age_min;
    ages_lim = age_max;
    min_no = 12;
    Xmin = 1;
    Xlim = ages_lim;
    
    viruses_ngs = viruses(find(viruses(:,ngs_index)>=7 & viruses(:,ngs_index)<=11 & viruses(:,ageid)<=ages_lim & viruses(:,ageid)>=ages_min),:);
    x_tot = viruses_ngs(:,index);
    x = x_tot;
    y_tot = viruses_ngs(:,Y_IDX);

    
    
    color_code = strvcat('G', 'B', 'R', 'K');
    for ngs_no=7:8
    ngs_no
    % Viral Strains with #NGS = 8
    dm1 = find(viruses_ngs(:,ngs_index)==ngs_no);
    dmI1 = zeros(size(x));
    dmI1(dm1)=1;
    
    viruses_ngs_no = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no & viruses_ngs(:,ageid)<=ages_lim),:);
    ac = group_aa_dist_byngs_byage(viruses_ngs_no);
    if isempty(ac)
    else
    siz = length(ac(:,1));
    for i=1:siz
      age = ac(i,1);
      charge = ac(i,3);
      amount = ac(i,4);
      if(amount<3)
        plot(age,charge,[color_code(ngs_no-6) 'o'], 'MarkerSize', 3, 'LineWidth', 1.2);
        %plot(age,charge,'Bo', 'MarkerSize', 2, 'MarkerFaceColor', 'B');
      elseif(amount<6)
        plot(age,charge,[color_code(ngs_no-6) 'o'], 'MarkerSize', 6, 'LineWidth', 1.2);
      else
          plot(age,charge,[color_code(ngs_no-6) 'o'], 'MarkerSize', 10);
      end
      hold on;
    end
    end
   

        
    if(length(viruses_ngs_no(:,1))>min_no)
        x1 = viruses_ngs_no(:,index);
        y = viruses_ngs_no(:,charge_index);
        
        % detect outliers
        reg_stats = regstats(y,x1);
        potential_outlier = reg_stats.cookd > (4/length(x1))*5;
        %scatter(x1(potential_outlier),y(potential_outlier), 'rx', 'LineWidth', 0.5)
        % drop outliers
        loc_x1 = find(viruses_ngs(:,ngs_index)==ngs_no);
        %x1(potential_outlier) = [];
        %y(potential_outlier)= [];
        viruses_ngs(loc_x1(potential_outlier),10) = 1;
        
        X = [ones(size(x1)) x1]; 
        [b,bint,r,rint,stats] = regress(y,X);
        x1fit = Xmin:ages_lim;
        YFIT = b(1) + b(2)*x1fit;
        plot(x1fit,YFIT, [color_code(ngs_no-6) '--'], 'LineWidth', 2);
        stats(2:3) 
        [slope pv se T] = getPfromTtest(X,y);
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
    end
    
    if exist('saveoutput')
        if strcmp(saveoutput,'s')
            save(['h1n1_ny_net_byage' num2str(ages_min) '-' num2str(ages_lim) '_reg_pvalue.mat'], 'pvalue', 'ptable', 'ages_min', 'ages_lim');
        end
    end
    
    pvalue = [];
    ptable = [];

    %create dummy variable
    viruses_ngs_noutlier = viruses_ngs(find(viruses_ngs(:,10)~=-1),:);
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
    Y1FIT = b(1) + b(end)*x1fit;
    plot(x1fit,Y1FIT, 'G-', 'LineWidth', 2);
    Y1FIT = b(1)+b(2) + b(end)*x1fit;
    plot(x1fit,Y1FIT, 'B-', 'LineWidth', 2);
    %Y1FIT = b(1)+b(3) + b(end)*x1fit;
    %plot(x1fit,Y1FIT, 'R-');
    %Y1FIT = b(1)+b(4) + b(end)*x1fit;
    %plot(x1fit,Y1FIT, 'K-');
    
    
    
    
    
    xlabel('host age(in yrs)');
    ylabel('net charge');
    xlim([0 round(Xlim./10)*10]);
    %set(gca,'XLim',[0 4])
    %set(gca,'XTick',[1 20 40 60])
    %set(gca,'XTickLabel',['1';' ';'20';' ';'40';' ';'3';' ';'4'])

    
    if exist('saveoutput')
        if strcmp(saveoutput,'s')
            save(['h1n1_ny_net_byage' num2str(age_min) '-' num2str(age_max) '_dummy_pvalue.mat'], 'pvalue', 'ptable', 'age_min', 'age_max');
        end
    end
%title([num2str(fix(ranges2(ageid,1))) '/' num2str(fix(ranges2(ageid,2)))]);
end
end
