% Plot the netcharge distribution by time
% Jul 23, 2013
% Hsiang-yu Yuan

n_tips=95; %total isolates
%proj = 'noram9506-20130730';
proj = 'nogap-20130731';
pairs = string2pairs(['dat/' proj '/hm_h1n1_noram_1995_2006.topo.nw.mcc.(time).trees'],n_tips);
load(['dat/' proj '/hm_h1n1_noram_charge_1995_2006_nogap.mat']);
load(['dat/' proj '/newick_elements.mat']);

%proj = 'noram9506-20130730';
%pairs = string2pairs(['dat/' proj '/hm_h1n1_noram_1995_2006.topo.nw.mcc.(time).trees'],n_tips);
%load(['dat/' proj '/hm_h1n1_noram_charge_1995_2006.mat']);
%load(['dat/' proj '/newick_elements.mat']);
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
dnum = datenum('2006/04/10','yyyy/mm/dd');
dpos = dnum./365.25;
rootdpos = dpos - distance_root(76);
distance_root_dpos = distance_root + rootdpos;

charge_name = [];
for id=1:length(charge_txt(:,1))
  loc = str2num(taxa(id).name);
  if isempty(loc) % This will be internal nodes
    loc = id; 
  end
  charge_name(loc,:) = charge_txt(id,:);
end

%charge_txt %read from hm_h3n2_ny_charge_1993.mat
nodes = 1:length(pairs);
[anc]= getParentList(pairs, 76); %Root for H1N1: 76 ABK79959
anc(find(anc==76)) = []; % remove latest strain from trunk lineage 
tr = nodes(ismember(nodes,anc));

hFig = figure;
Figw = 1080;
Figh = 260;
set(hFig, 'Position', [100 100 Figw Figh]);
for ngs_no=7:8
ngs_no
% Total nodes with NGS=8
ngs = find(charge_name(:,5)==ngs_no); 
% find the trunk nodes
ngs_tr = tr(ismember(tr,ngs));  
% find the non-trunk nodes
ngs_nt = ngs(~ismember(ngs, anc));
% find the external tips
ngs_e = ngs(ismember(ngs, [1:n_tips]));
% find the internal tips of non-trunk
ngs_i = intersect([n_tips+1:n_total], ngs_nt);


% plot
subplot(2,1,ngs_no-6);
hold on;
plot(distance_root_dpos(ngs_tr), charge_name(ngs_tr,1),'ro');
plot(distance_root_dpos(ngs_i), charge_name(ngs_i,1),'b.');
plot(distance_root_dpos(ngs_e), charge_name(ngs_e,1),'x', 'Color', [0 .8 0]);
if (ngs_no == 7)
    line([1996.2 1996.2],[2 10], 'Color', [.8 .8 .8]);
    %line([2001.6 2001.6],[2 10], 'Color', [.8 .8 .8]);
end

xlim([rootdpos rootdpos+distance_root(76)]);
if ngs_no < 8
    set(gca,'XTickLabel',['']) % remove tick label
end
ylabel('net charge');
ylim([3.5 8.5]);
end

xlabel('year of isolations');


