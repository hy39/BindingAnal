function [ netcharge ] = calNetCharge( )
n_tips = 686;
n_total = n_tips*2-1;

%pairs = string2pairs('dat/ancestor_20130702/hm_h3n2_ny_dna_beast_1993_2006.topo.nw.mcc.(time).trees',665);
pairs = string2pairs('dat/ancestor_20130713/hm_h3n2_ny_dna_beast_1993.topo.nw.mcc.(time).trees',686);
%[anc]= getParentList(pairs, 618); %618 ABJ53482_CY016995_2006/04/06_11_17,
%anc(find(anc==618)) = [];
[anc]= getParentList(pairs, 633); %633 ABJ53482_CY016995_2006/04/06_11_17,
anc(find(anc==633)) = [];

%load('dat/ancestor_20130702/hm_h3n2_ny_charge_1993_2006.mat');
load('dat/ancestor_20130713/hm_h3n2_ny_charge_1993.mat');
load('dat/ancestor_20130713/newick_elements.mat'); % Read taxa

charge_name = []; % Use original lable name not the sequential ID, although it give same results. 
for id=1:length(charge_txt(:,1))
  loc = str2num(taxa(id).name);
  if isempty(loc) % This will be internal nodes
    loc = id; 
  end
  charge_name(loc,:) = charge_txt(id,:);
end

% find the non trunk
nodes = [1:n_total];
non_trunk = nodes(~ismember(1:n_total, anc))

ngs8 = find(charge_name(:,5)==8);
ngs9 = find(charge_name(:,5)==9);
ngs10 = find(charge_name(:,5)==10);
ngs11 = find(charge_name(:,5)==11);

% find the non-trunk nodes
ngs8nt = ngs8(~ismember(ngs8, anc));
ngs9nt = ngs9(~ismember(ngs9, anc));
ngs10nt = ngs10(~ismember(ngs10, anc));
ngs11nt = ngs11(~ismember(ngs11, anc));
ngsnt = {ngs8nt, ngs9nt, ngs10nt, ngs11nt};

% find the external tips
ngs8e = ngs8(ismember(ngs8, [1:n_tips]));
ngs9e = ngs9(ismember(ngs9, [1:n_tips]));
ngs10e = ngs10(ismember(ngs10, [1:n_tips]));
ngs11e = ngs11(ismember(ngs11, [1:n_tips]));
ngse = {ngs8e, ngs9e, ngs10e, ngs11e};

% find the internal tips of non-trunk
ngs8i = intersect([n_tips+1:n_total], ngs8nt);
ngs9i = intersect([n_tips+1:n_total], ngs9nt);
ngs10i = intersect([n_tips+1:n_total], ngs10nt);
ngs11i = intersect([n_tips+1:n_total], ngs11nt);
ngsi = {ngs8i, ngs9i, ngs10i, ngs11i};

% find the trunk nodes
anc8 = anc(ismember(anc,ngs8));
anc9 = anc(ismember(anc,ngs9));
anc10 = anc(ismember(anc,ngs10));
anc11 = anc(ismember(anc,ngs11));
ngstk = {anc8, anc9, anc10, anc11};

%netcharge = {[]};
cnt = 8;
for n=1:4
netcharge(n).ngs = cnt;
netcharge(n).trunk_mean = mean(charge_name(cell2mat(ngstk(n)),1));
netcharge(n).trunk_var = var(charge_name(cell2mat(ngstk(n)),1));
netcharge(n).trunk_sd = (var(charge_name(cell2mat(ngstk(n)),1))).^0.5;
netcharge(n).nontrunk_mean = mean(charge_name(cell2mat(ngsnt(n)),1));
netcharge(n).nontrunk_var = var(charge_name(cell2mat(ngsnt(n)),1));
netcharge(n).nontrunk_sd = (var(charge_name(cell2mat(ngsnt(n)),1))).^0.5;
netcharge(n).internal_mean = mean(charge_name(cell2mat(ngsi(n)),1));
netcharge(n).internal_var = var(charge_name(cell2mat(ngsi(n)),1));
netcharge(n).internal_sd = (var(charge_name(cell2mat(ngsi(n)),1))).^0.5;
netcharge(n).external_mean = mean(charge_name(cell2mat(ngse(n)),1));
netcharge(n).external_var = var(charge_name(cell2mat(ngse(n)),1));
netcharge(n).external_sd = (var(charge_name(cell2mat(ngse(n)),1))).^0.5;
cnt = cnt + 1;
end



%mean(charge_txt(anc9,1))
%var(charge_txt(anc9,1))
%mean(charge_txt(ngs9a,1))
%var(charge_txt(ngs9a,1))

%mean(charge_txt(anc10,1))
%var(charge_txt(anc10,1))
%mean(charge_txt(ngs10a,1))
%var(charge_txt(ngs10a,1))

%mean(charge_txt(anc11,1))
%var(charge_txt(anc11,1))
%mean(charge_txt(ngs11a,1))
%var(charge_txt(ngs11a,1))




