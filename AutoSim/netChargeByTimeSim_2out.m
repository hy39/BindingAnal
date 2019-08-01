
% NOTE: this is adapted from plotNetchargeByTimeSim.m

function[intMat, extMat] = netChargeByTimeSim_2out(infileV, infileH, intMat, extMat)

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

% for each node, find its ancestor
charge_name = [charge_name zeros(size(charge_name,1),2)];
for id = 1:size(charge_name,1)
    % 1. find the current node's row position in the table b <-> its parent
    ancestor = find(b(:,1:2)==id);
    % 2. enter the parent's charge on the 2nd column
    if (isempty(ancestor))
        charge_name(id,2) = NaN;
    else
        charge_name(id,2) = charge_name(ancestor,1);
    end
    % 3. check whether the current node belongs to ngs_i <-> internal
    % node, and put a marker on the 3rd column accordingly
    if (ismember(id,ngs_i))
        charge_name(id,3) = 'i';
    else charge_name(id,3) = 'e';
    end
end

% 4. find the average change in charge for internal and external nodes and
% append the values to intMat and extMat respectively. NOTE: for now,
% intMat and extMat are designed to each have a dimension of 200x599 (200
% simulations x 300+299 external & internal nodes respectively per
% simulation, with many entries left as NaN in each matrix).
int_count = 0; ext_count = 0;
tempIntMat = zeros(1,n_tips); tempExtMat = zeros(1,n_tips);
for i = 1:size(charge_name,1)
    if(charge_name(i,3)=='i')
        tempIntMat(int_count+1) = charge_name(i,1) - charge_name(i,2);
        int_count=int_count+1;
    else
        tempExtMat(ext_count+1) = charge_name(i,1) - charge_name(i,2);
        ext_count=ext_count+1;
    end
end
emptyRows = n_tips - int_count; % n_tips == smpno
if (emptyRows>0)
    tempIntMat(int_count:end) = NaN;
end
emptyRows = n_tips - ext_count;
if (emptyRows>0)
    tempExtMat(ext_count:end) = NaN;
end

intMat = [intMat; tempIntMat];
extMat = [extMat; tempExtMat];

end


