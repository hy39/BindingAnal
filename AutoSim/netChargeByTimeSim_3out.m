
% NOTE: this is adapted from plotNetchargeByTimeSim.m

function[VintMat, VextMat, VoptMat] = netChargeByTimeSim_3out(infileV, infileH, VintMat, VextMat, VoptMat)

disp(['loading file ' infileV]);
load(infileV);
% import b, d, names
n_tips=length(b(:,1))+1;
b(:,3)=n_tips+1:n_tips*2-1;
distance_root = [];
pairs = zeros(n_tips*2-1,1);
for i=1:n_tips-1
  pairs(b(i,1))=b(i,3);
  pairs(b(i,2))=b(i,3);
end
n_total = length(pairs);
disp(['total number of tips: ' num2str(n_tips)]);
disp(['total number of tips and internal nodes: ' num2str(n_total)]); 

% find the all descendents for a given node anc_node
des_list = [];
%anc_node = 599; %20160122
anc_node = n_tips*2-1;
des_list1 = anc_node;
i = 1;
while i<length(des_list1)+1 
    des_node = find(pairs(:,1) == des_list1(i));
    des_list1 = [des_list1; des_node];
    i = i + 1;
end
pairs(des_list1,2)=1; % the first clade

% calculate the distance to root for each node
disp(['calculating distance to root..']);
for i=1:length(pairs)
 [lin]= getParentList(pairs, i);
 distance = 0;
 for j=1:length(lin)
    if lin(j) == 0
        break;
    end
    lin(j);
    distance = distance + d(lin(j)); 
 end
 distance_root(i,:) = distance;
end
root_start = 0;
distance_root = distance_root + root_start;

% extract binding scores. 
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

% find the 1st main trunk
nodes = 1:length(pairs);
[anc1]= getParentList(pairs, n_tips); % the id of the latest strain in the tree is n_tips 
anc1(find(anc1==n_tips)) = []; % remove latest strain from trunk lineage 
tr1 = nodes(ismember(nodes,anc1));

% if needed, find the 2nd main trunk 
%[anc2]= getParentList(pairs, 291);
%anc2(find(anc2==291)) = []; % remove latest strain from trunk lineage 
%tr2 = nodes(ismember(nodes,anc2));

ngs_total = 1:n_total;
Lia = ismember(ngs_total,des_list1);
ngs = ngs_total(Lia);
% find the trunk nodes
ngs_tr = tr1(ismember(tr1,ngs));  
% find the non-trunk nodes
ngs_nt = ngs(~ismember(ngs, anc1));
% find the external tips
ngs_e = ngs(ismember(ngs, [1:n_tips]));
% find the internal tips
ngs_i = [n_tips+1:n_total];
% find the internal tips of non-trunk
ngs_int = intersect([n_tips+1:n_total], ngs_nt);

%plot internal nodes V and external tips V
set_dpos_in = distance_root_dpos(ngs_i);
set_charge_in = charge_name(ngs_i,1);
set_dpos_ext = distance_root_dpos(ngs_e);
set_charge_ext = charge_name(ngs_e,1);
Vint = 0;
Vext = 0;
interval = VintMat(2,1) - VintMat(1,1);
for i=1:interval:round(max(set_dpos_in))
%for i=1:30:round(max(set_dpos_in)) % use 1:15 as previous setting
  Vint = [Vint mean(set_charge_in(find(set_dpos_in>i & set_dpos_in<i+2*interval)))]; % NOTE: the 2*interval is to smoothen the curve
end
if (length(Vint) > size(VintMat,1))
    Vint = Vint(1:size(VintMat,1));
elseif (length(Vint) < size(VintMat,1))
    Vint = [Vint zeros(1, size(VintMat,1)-length(Vint))];
end
%Vint(1) = [];
VintMat = [VintMat Vint'];
for i=1:interval:round(max(set_dpos_ext))
%for i=1:30:round(max(set_dpos_ext)) % use 1:30 as previous setting
  Vext = [Vext mean(set_charge_ext(find(set_dpos_ext>i & set_dpos_ext<i+2*interval)))];
end
if (length(Vext) > size(VextMat,1))
    Vext = Vext(1:size(VextMat,1));
elseif (length(Vext) < size(VextMat,1))
    Vext = [Vext zeros(1, size(VextMat,1)-length(Vext))];
end
%Vext(1) = [];
VextMat = [VextMat Vext'];
%%%shift=10;
%%%plot(1+shift:15:length(Vint)*15+shift, Vint, '-');
%%%plot(1+shift:15:length(Vext)*15+shift, Vext, '-');

% Voptimum
load(infileH);
totalDays = VoptMat(end,1);
interval = VoptMat(2,1)-VoptMat(1,1);
shift = VoptMat(1,1);
V = Vopt(1:totalDays/interval)';
if (length(V)>size(VoptMat,1))
    V = V(1:size(VoptMat,1));
end
VoptMat = [VoptMat V];
%%%plot(1+shift:5:totalDays+shift,Vopt(1:totalDays/5),'-','Color', [0.4 0.4 0.4]);

end
