% Plot the netcharge distribution and calculate diversity by time
% May 29, 2014
%
% Extract the mutational effect from parent nodes
% May 22, 2013
% Hsiang-yu Yuan

n_tips=686; %total isolates
pairs = string2pairs('dat/ancestor_20130713/hm_h3n2_ny_dna_beast_1993.topo.nw.mcc.(time).trees',n_tips);
load('dat/ancestor_20130713/hm_h3n2_ny_charge_1993_2006.mat');
load('dat/ancestor_20130713/newick_elements.mat');
distance_root = [];
n_total = length(pairs);

for i=1:length(pairs)
 [lin]= getParentList(pairs, i);
distance = 0;
for j=1:length(lin)
if lin(j) == 0
  break;
end
distance = distance + b(lin(j)); 
end
distance_root(i,:) = distance;
end

root_start = 0.15;
distance_root = distance_root + root_start;

% update distance_to_root to dates
dnum = datenum('2006/04/06','yyyy/mm/dd');
dpos = dnum./365.25;
rootdpos = dpos - distance_root(633);
distance_root_dpos = distance_root + rootdpos;

charge_name = [];
for id=1:length(charge_txt(:,1))
  loc = str2num(taxa(id).name);
  if isempty(loc) % This will be internal nodes
    loc = id; 
  end
  charge_name(loc,:) = charge_txt(id,:);
end

%%% 20140522
%%% Add mutational effect to charge_name array
charge_name = [charge_name zeros(n_total,1)];
for loc=1:length(charge_name(:,1))
  parent = pairs(loc);
  if parent~=0
   mutation = charge_name(loc,1) - charge_name(parent,1);
   charge_name(loc,6) = mutation;
  end
end

%charge_txt %read from hm_h3n2_ny_charge_1993.mat
nodes = 1:length(pairs);
[anc]= getParentList(pairs, 633);
anc(find(anc==633)) = []; % remove latest strain from trunk lineage 
tr = nodes(ismember(nodes,anc));


ngs_no=[8:11];
ngs_no
ngs_idx = [];
% Total nodes with NGS=8
for ng = ngs_no(1):ngs_no(end)
  
    ngs_tmp = find(charge_name(:,5)==ng);
    ngs_idx = [ngs_idx; ngs_tmp]; 
end
% find the trunk nodes
ngs_tr = tr(ismember(tr,ngs_idx));  
% find the non-trunk nodes
ngs_nt = ngs_idx(~ismember(ngs_idx, anc));
% find the external tips
ngs_e = ngs_idx(ismember(ngs_idx, [1:n_tips]));
% find the internal tips of non-trunk
ngs_i = intersect([n_tips+1:n_total], ngs_nt);
%ngs_ia = [ngs_i;ngs_tr'];
ngs_ia = [n_tips+1:n_total];
mean(abs(charge_name(ngs_ia,6)))
mean(abs(charge_name(ngs_e,6)))





