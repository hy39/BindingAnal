function [] = plot_net_byage( age_min, age_max )

p = path
path(p,'..');
p = path
path(p,'../lib');
p = path
path(p,'../../lib');
addpath('Z:\Projects\Binding_Evo\NetCharge\mat\netcharge2')
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
ranges2 = [1993.601 2005.6];   
%load('../dat/h3n2_ny/hm_h3n2_ny_merged_data.mat');
load('../../dat/h3n2_ny/hm_h3n2_ny_merged_data.mat');
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

    figure;
    hold;
    
    ages_min = age_min;
    ages_lim = age_max;
    min_no = 12;
    
    viruses_ngs = viruses(find(viruses(:,ngs_index)>=8 & viruses(:,ngs_index)<=11 & viruses(:,ageid)<=ages_lim & viruses(:,ageid)>=ages_min),:);
    x_tot = viruses_ngs(:,index);
    y_tot = viruses_ngs(:,Y_IDX);
    x = x_tot;
    
    
    color_code = strvcat('G', 'B', 'C', 'K');
    for ngs_no=8:11
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
        plot(age,charge,[color_code(ngs_no-8) 'o'], 'MarkerSize', 3);
        %plot(age,charge,'Bo', 'MarkerSize', 2, 'MarkerFaceColor', 'B');
      elseif(amount<8)
        plot(age,charge,[color_code(ngs_no-8) 'o'], 'MarkerSize', 8);
      else
          plot(age,charge,[color_code(ngs_no-8) 'o'], 'MarkerSize', 12);
      end
      hold on;
    end
    end
    
    if(length(viruses_ngs_no(:,1))>min_no)
    x = viruses_ngs_no(:,index);
    y = viruses_ngs_no(:,charge_index);
        X = [ones(size(x)) x]; 
        [b,bint,r,rint,stats] = regress(y,X);
        x1fit = 0:ages_lim;
        YFIT = b(1) + b(2)*x1fit;
        plot(x1fit,YFIT, [color_code(ngs_no-8) '-']);
        stats(2:3) 
        [slope pv se T] = getPfromTtest(X,y);
        pv = round(pv.*100)./100;
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
    %ptable(j1,:) = [pvalue(j1).season pvalue(j1).ngs pvalue(j1).slope bint(2,1) bint(2,2) pvalue(j1).stats(1) pvalue(j1).stats(2) pvalue(j1).stats(3)];
    ptable(j1,:) = [pvalue(j1).season pvalue(j1).ngs slope T pv se stats(1)];
    end
    end
    
    
    % Viral Strains with #NGS = 9
    ngs_no = 9;
    dm2 = find(viruses_ngs(:,ngs_index)==ngs_no);
    dmI2 = zeros(size(x));
    dmI2(dm2)=1;
    
    viruses_ngs9 = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no & viruses_ngs(:,ageid)<=ages_lim),:);
    ac = group_aa_dist_byngs_byage(viruses_ngs9);
    if isempty(ac)
    else
    siz = length(ac(:,1));
    for i=1:siz
      age = ac(i,1);
      charge = ac(i,3);
      amount = ac(i,4);
      if(amount<3)
        plot(age,charge,'Bo', 'MarkerSize', 5);
        %plot(age,charge,'Bo', 'MarkerSize', 2, 'MarkerFaceColor', 'B');
      elseif(amount<8)
        plot(age,charge,'Bo', 'MarkerSize', 8);
      else
          plot(age,charge,'Bo', 'MarkerSize', 12);
      end
      hold on;
    end
    end
     
    if(length(viruses_ngs9(:,1))>min_no)
        x = viruses_ngs9(:,index);
        y = viruses_ngs9(:,charge_index);
        X = [ones(size(x)) x]; 
        [b,bint,r,rint,stats] = regress(y,X);
        x1fit = 0:ages_lim;
        YFIT = b(1) + b(2)*x1fit;
        plot(x1fit,YFIT, 'B-');
        stats(2:3) 
        [slope pv se T] = getPfromTtest(X,y);
        pv = round(pv.*100)./100;
        slope = round(slope.*1000.*100)./100;
        se = round(se.*1000.*100)./100;
    j1 = length(pvalue)+1;
    %pvalue(j1).metadata = 'season start, season end, slope, upper bound, lower bound, R2, T score, p value, standard error';
    pvalue(j1).metadata = 'season start, season end, slope, T score, p value, standard error, R2';
    pvalue(j1).season = ranges2;
    pvalue(j1).ngs = ngs_no;
    pvalue(j1).slope = b(2);
    pvalue(j1).bint = bint;
    pvalue(j1).stats = stats;
    pvalue(j1).p = stats(3);
    ptable(j1,:) = [pvalue(j1).season pvalue(j1).ngs slope T pv se stats(1)];
    end
    
    % Viral Strains with #NGS = 10
    ngs_no = 10;
    dm3 = find(viruses_ngs(:,ngs_index)==ngs_no);
    dmI3 = zeros(size(x));
    dmI3(dm3)=1;
    
    viruses_ngs10 = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no & viruses_ngs(:,ageid)<=ages_lim),:);
    ac = group_aa_dist_byngs_byage(viruses_ngs10);
    if isempty(ac)
    else
    siz = length(ac(:,1));
    for i=1:siz
      age = ac(i,1);
      charge = ac(i,3);
      amount = ac(i,4);
      if(amount<3)
        plot(age,charge,'Co', 'MarkerSize', 5);
        %plot(age,charge,'Bo', 'MarkerSize', 2, 'MarkerFaceColor', 'B');
      elseif(amount<8)
        plot(age,charge,'Co', 'MarkerSize', 8);
      else
          plot(age,charge,'Co', 'MarkerSize', 12);
      end
      hold on;
    end
    end
    
    if(length(viruses_ngs10(:,1))>min_no)
        x = viruses_ngs10(:,index);
        y = viruses_ngs10(:,charge_index);
        X = [ones(size(x)) x]; 
        [b,bint,r,rint,stats] = regress(y,X);
        x1fit = 0:ages_lim;
        YFIT = b(1) + b(2)*x1fit;
        plot(x1fit,YFIT, 'C-');
        stats(2:3) 
        [slope pv se T] = getPfromTtest(X,y);
        pv = round(pv.*100)./100;
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
    ptable(j1,:) = [pvalue(j1).season pvalue(j1).ngs slope T pv se stats(1)];  
    end
    
    % Viral Strains with #NGS = 11
    ngs_no = 11;
    dm4 = find(viruses_ngs(:,ngs_index)==ngs_no);
    dmI4 = zeros(size(x));
    dmI4(dm4)=1;
    viruses_ngs11 = viruses_ngs(find(viruses_ngs(:,ngs_index)==ngs_no & viruses_ngs(:,ageid)<=ages_lim),:);
    ac = group_aa_dist_byngs_byage(viruses_ngs11);
    if isempty(ac)
    else
    siz = length(ac(:,1));
    for i=1:siz
      age = ac(i,1);
      charge = ac(i,3);
      amount = ac(i,4);
      if(amount<3)
        plot(age,charge,'Ko', 'MarkerSize', 5);
        %plot(age,charge,'Bo', 'MarkerSize', 2, 'MarkerFaceColor', 'B');
      elseif(amount<8)
        plot(age,charge,'Ko', 'MarkerSize', 8);
      else
          plot(age,charge,'Ko', 'MarkerSize', 12);
      end
      hold on;
    end
    end
    
    if(length(viruses_ngs11(:,1))>min_no)
        x = viruses_ngs11(:,index);
        y = viruses_ngs11(:,charge_index);
        X = [ones(size(x)) x]; 
        [b,bint,r,rint,stats] = regress(y,X);
        x1fit = 0:ages_lim;
        YFIT = b(1) + b(2)*x1fit;
        plot(x1fit,YFIT, 'K-');
        stats(2:3) 
        [slope pv se T] = getPfromTtest(X,y);
        pv = round(pv.*100)./100;
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
    ptable(j1,:) = [pvalue(j1).season pvalue(j1).ngs slope T pv se stats(1)];
    end 
    
    save('h3n2_ny_netcharge_ages_reg_pvalue.mat', 'pvalue', 'ptable');


%axis([ 0, 0, 14, 20])

%set(gca,'XLim',[0 ages_lim])
%set(gca,'YLim',[12 21]) %set the range
%pos1 = [100   100   1480   762];
%set(hfig,'Position',pos1); 


%set(gca,'YLim',[68 75]) %set the range
%x_scale = [1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006];
%x_label = ['1995/1996'; '1996/1997'; '1997/1998'; '1998/1999'; '1999/2000'; '2000/2001'; '2001/2002'; '2002/2003'; '2003/2004'; '2004/2005'; '2005/2006'];
%set(gca,'XTick',x_scale);
%set(gca,'XTickLabel',x_label);
title([num2str(fix(ranges2(ageid,1))) '/' num2str(fix(ranges2(ageid,2)))]);
end
end
