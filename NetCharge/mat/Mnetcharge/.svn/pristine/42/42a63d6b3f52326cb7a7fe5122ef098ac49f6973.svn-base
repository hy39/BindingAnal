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
hfig = figure;
i = 1;


y_lim = [30 45];

pvalue = [];
ptable = [];


    ngs_index = 4;
    charge_index = 3; %netcharge
    %charge_index = 5; %total charged amino acids
    ageid = 1;
    kid = 6;
    index = charge_index; % 1:age, 2:isolation date, 6:times of previouis infection, 3:netcharge 
    bding_charge_index = 7;
    bding_total_index = 9;
    Y_IDX = bding_total_index;
    min_no = 12;
    
    % #NGS=8
    TF = find(Viruses(:,2)>ranges2(1,1) & Viruses(:,2)<ranges2(1,2) &  Viruses(:,1)~=0);
    viruses = Viruses(TF,:);

    
    dm_ngs = dominance(i,:);
    charge_lim = 20;
    ages_lim = 70;
    min_no = 12;
    Xlim = charge_lim;
    Xmin = 14;
    
    
    
    for ngs_no=8:11
        
    if ngs_no==8
        subplot(2,2,1);
        set(gca,'YLim',y_lim(i,:));
        hold;
    elseif ngs_no==9
        subplot(2,2,2);
        set(gca,'YLim',y_lim(i,:));
        hold;
    elseif ngs_no==10
        subplot(2,2,3);
        set(gca,'YLim',y_lim(i,:));
        hold;
    elseif ngs_no==11
        subplot(2,2,4);
        set(gca,'YLim',y_lim(i,:));
        hold;
    end
    set(gca,'XLim',[14 20])
    title(['#NGS=' num2str(ngs_no)]);
    
    %viruses_ngs = viruses(find(viruses(:,ngs_index)>=8 & viruses(:,ngs_index)<=11 & viruses(:,ageid)<=ages_lim),:);
    viruses_ngs = viruses(find(viruses(:,ngs_index)==ngs_no & viruses(:,ageid)<=ages_lim),:);
    x = viruses_ngs(:,index);
    y = viruses_ngs(:,Y_IDX);
    %ngs_no = 8;
    %dm1 = find(viruses_ngs(:,ngs_index)==ngs_no);
    %dmI1 = zeros(size(x));
    %dmI1(dm1)=1;
    
    %ngs_no = 9;
    %dm2 = find(viruses_ngs(:,ngs_index)==ngs_no);
    %dmI2 = zeros(size(x));
    %dmI2(dm2)=1;

    %ngs_no = 10;
    %dm3 = find(viruses_ngs(:,ngs_index)==ngs_no);
    %dmI3 = zeros(size(x));
    %dmI3(dm3)=1;
    
    %ngs_no = 11;
    %dm4 = find(viruses_ngs(:,ngs_index)==ngs_no);
    %dmI4 = zeros(size(x));
    %dmI4(dm4)=1;
    
    %X = [ones(size(x)) dmI1 dmI2 dmI3 x];
    X = [ones(size(x)) x]; 
    [b,bint,r,rint,stats] = regress(y,X);
    [slope pv se T] = getPfromTtest(X,y);
       
    
    j1 = length(pvalue)+1;
    pvalue(j1).metadata = 'season start, season end, slope, T score, p value, standard error, R2';
    pvalue(j1).season = ranges2;
    pvalue(j1).ngs = ngs_no;
    pvalue(j1).slope = slope;
    ptable(j1,:) = [pvalue(j1).season ngs_no slope T pv se stats(1)];
    
    if ngs_no == 8
        color = 'G.';
    elseif ngs_no == 9
        color = 'B.';
    elseif ngs_no == 10
        color = 'C.';
    elseif ngs_no == 11
        color = 'K.';
    end
        plot(x, y, color);
        x1fit = Xmin:Xlim;
        YFIT = b(1) + b(2)*x1fit;
        plot(x1fit,YFIT, 'color', [0.5 0.5 0.5]);
    
    end
    
    save('h3n2_ny_bding_netcharge_allyears_pvalue.mat', 'pvalue', 'ptable');
    
  

