% Plot the netcharge distribution by time
% Jul 23, 2013
% Hsiang-yu Yuan

n_tips=686; %total isolates
pairs = string2pairs('dat/ancestor_20130713/hm_h3n2_ny_dna_beast_1993.topo.nw.mcc.(time).trees',n_tips);
load('dat/ancestor_20130713/hm_h3n2_ny_charge_1993_2006.mat');
load('dat/ancestor_20130713/hm_h3n2_ny_binding_1993.mat');
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
  charge_name(loc,1:5) = charge_txt(id,1:5);
  charge_name(loc,6) = binding_txt(id,2);
end


%%% 20140522
%%% Add mutational effect to charge_name array
charge_name = [charge_name zeros(n_total,1)];
for loc=1:length(charge_name(:,1))
  parent = pairs(loc);
  if parent~=0
   mutation = charge_name(loc,6) - charge_name(parent,6);
   charge_name(loc,7) = mutation;
  end
end


%charge_txt %read from hm_h3n2_ny_charge_1993.mat
nodes = 1:length(pairs);
[anc]= getParentList(pairs, 633);
anc(find(anc==633)) = []; % remove latest strain from trunk lineage 
tr = nodes(ismember(nodes,anc));

hFig = figure;
Figw = 1080;
Figh = 400;
set(hFig, 'Position', [100 100 Figw Figh]);
for ngs_no=8:11
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
subplot(4,1,ngs_no-7);
hold on;
% Both trunk and internal nodes use circle, external nodes use x.
plot(distance_root_dpos(ngs_tr), charge_name(ngs_tr,6),'ro');
plot(distance_root_dpos(ngs_i), charge_name(ngs_i,6),'b.');
plot(distance_root_dpos(ngs_e), charge_name(ngs_e,6),'x', 'Color', [0 .8 0]);
%xlim([0 13]);
if (ngs_no == 8)
    line([1998.2 1998.2],[25 50], 'Color', [.5 .5 .5]);
end
if (ngs_no == 9)
    line([1995.2 1995.2],[25 45], 'Color', [.5 .5 .5]);
end
if (ngs_no == 10)
    line([2000.18 2000.18],[25 45], 'Color', [.5 .5 .5]);
    line([2003.25 2003.25],[25 45], 'Color', [.5 .5 .5]);
    line([2004.05 2004.05],[25 45], 'Color', [.5 .5 .5]);
end
if (ngs_no == 11)
    line([2000.05 2000.05],[25 45], 'Color', [.5 .5 .5]);
    line([2004.1 2004.1],[25 45], 'Color', [.5 .5 .5]);
end

xlim([rootdpos rootdpos+distance_root(633)]);
if ngs_no < 11
    set(gca,'XTickLabel',['']) % remove tick label
end

ylabel('binding score');
ylim([29 44]);
end

xlabel('time to the root (years)');


