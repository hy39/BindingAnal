% Plot the netcharge distribution by time
% Phylogeny: indiv_genealogy_300.nx.branch.tree
% Feb 10, 2014
% Hsiang-yu Yuan


n_tips=300; %total isolates
%pairs = string2pairs('dat/ancestor_20130713/hm_h3n2_ny_dna_beast_1993.topo.nw.mcc.(time).trees',n_tips);
%%pairs = string2pairs('dat/simphylo/genealogy_300(2).nw',n_tips);
%load('dat/ancestor_20130713/hm_h3n2_ny_charge_1993_2006.mat');
%load('dat/ancestor_20130713/newick_elements.mat');
%load('dat/simphylo/differential/tree_300.mat');
load('dat/simphylo/equal/tree_300.mat');
b(:,3)=n_tips+1:n_tips*2-1;
distance_root = [];
pairs = zeros(n_tips*2-1,1);
for i=1:n_tips-1
  pairs(b(i,1))=b(i,3);
  pairs(b(i,2))=b(i,3);
end
n_total = length(pairs);

    
%%
%% 09 Feb 2014
%% I need a subfunction that can find all the subset from an internal node.
%% Then I can define a subset of the viral clade.
%% Hint: use pairs. E.g., if the ancestor node is 570, then use find(pairs == 570) to extract the 
%% decsendants.
des_list = [];
anc_node = 599;
des_list1 = anc_node;
i = 1;
while i<length(des_list1)+1 
    des_node = find(pairs(:,1) == des_list1(i));
    des_list1 = [des_list1; des_node];
    i = i + 1;
end
pairs(des_list1,2)=1; % the first clade


for i=1:length(pairs)
 [lin]= getParentList(pairs, i);
 distance = 0;
 for j=1:length(lin)
    if lin(j) == 0
        break;
    end
    lin(j)
    distance = distance + d(lin(j)); 
 end
 distance_root(i,:) = distance;
end

root_start = 0;
distance_root = distance_root + root_start;

%%%
%%%
% 07 Feb 2014, Need to use names object to extract binding scores. 
% update distance_to_root to dates
dnum = datenum('2006/04/06','yyyy/mm/dd');
dpos = dnum./365.25;
%rootdpos = dpos - distance_root(599);
rootpos = 0;
distance_root_dpos = distance_root;

charge_name = [];
for id=1:length(names(1,:))
   str = names(1,id);
   label = '\d+[';
   x1 = regexp(str, label, 'match');
   x1_str = char(x1{1}); 
   x1_str = x1_str(1:end-1);
   x1 = str2num(x1_str);
  
   %%%vvvvv still working for this
   bind = 'Bindingscr.set="(\d+.\d+)"';
   x2 = regexp(str, bind, 'match');
   x2_str = char(x2{1});
   bindingscore = regexp(x2_str,'\d+.\d+','match');
   bindingscore = str2num(char(bindingscore{1}));
%  if isempty(loc) % This will be internal nodes
%    loc = id; 
%  end
   charge_name(id,1) = bindingscore;
end

% set the latest strain for each clade 
nodes = 1:length(pairs);
[anc1]= getParentList(pairs, 300);
anc1(find(anc1==300)) = []; % remove latest strain from trunk lineage 
tr1 = nodes(ismember(nodes,anc1));

%[anc2]= getParentList(pairs, 291);
%anc2(find(anc2==291)) = []; % remove latest strain from trunk lineage 
%tr2 = nodes(ismember(nodes,anc2));



hFig = figure;
Figw = 1080;
Figh = 400;
set(hFig, 'Position', [100 100 Figw Figh]);

%for clade 1 to 3
%ngs_no
% Total nodes with clade=1

ngs_total = 1:n_total;
Lia = ismember(ngs_total,des_list1);
ngs = ngs_total(Lia);
% find the trunk nodes
ngs_tr = tr1(ismember(tr1,ngs));  
% find the non-trunk nodes
ngs_nt = ngs(~ismember(ngs, anc1));
% find the external tips
ngs_e = ngs(ismember(ngs, [1:n_tips]));
% find the internal tips of non-trunk
ngs_i = intersect([n_tips+1:n_total], ngs_nt);

%% calculate the average binding avidity change
charge_diff = 0;
charge_diff_e = 0;
for i=1:length(ngs_i)
  charge_diff(i) = abs(charge_name(ngs_i(i))-charge_name(pairs(ngs_i(i))));
end

for e=1:length(ngs_e)
  charge_diff_e(e) = abs(charge_name(ngs_e(e))-charge_name(pairs(ngs_e(e))));
end

% plot
%subplot(2,1,1);
plot(distance_root_dpos(ngs_tr), charge_name(ngs_tr,1),'b.');
hold on;
plot(distance_root_dpos(ngs_i), charge_name(ngs_i,1),'b.');
plot(distance_root_dpos(ngs_e), charge_name(ngs_e,1),'x', 'Color', [0 .8 0]);


ylabel('binding score');
xlabel('time to the root');


