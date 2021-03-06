function void = BuildTree_indiv_nexus(filename_infectionTreeData, n_seqs)
% Produce the tree file as nexus format
p = path;
p = path(p,'lib/');
global epi_params;

%load(filename_simData); % Is this line needed? 20130716
load(filename_infectionTreeData);
who

folder = char(regexp(filename_infectionTreeData,'.+/','match'));
pos_s = regexp(filename_infectionTreeData,'/');
if isempty(pos_s)
    pos_s = 0;
end
filename = filename_infectionTreeData(pos_s(end)+1:end);
outfile_treeData = strcat(folder, 'indiv_treeData_', int2str(n_seqs));            
outfile_tree = strcat(folder, filename, '.tree');
n_tot_samples = n_seqs; 

% seq_times are the death times of each individuals
% indiv_sampled stored the index (the ith) of the sampled viruses
[seq_times, indiv_sampled] = GetIndividualsSampled_indiv(epi_params, births, deaths, n_tot_samples);

% Should retrieve viruses binding avidities Jul 17, 2013
% Access virus_traits.mat
indiv_sampled_binding_s = binding(indiv_sampled);
indiv_sampled_infectionK = infectionK(indiv_sampled);
if strcmp(epi_params.lg,'c')
    indiv_sampled_antigenicDrift = antigenicTot(indiv_sampled);
end
if strcmp(epi_params.lg,'m')
    indiv_sampled_antigenicDrift = zeros(length(indiv_sampled),1);
end

indiv_sampled_infecteddays = seq_times' - births(indiv_sampled);
indiv_sampled_binding = [];
annotation = epi_params.annotation;

% skip this part. only retrive initial virus binding avidity 
%for i = 1:n_tot_samples
  %indiv_sampled_binding(i) = 0.5;
  %indiv_sampled_binding(i) = getVChange_ode(indiv_sampled_binding_s(i),indiv_sampled_infectionK(i)-1,indiv_sampled_infecteddays(i)); %getVChange_ode(initialV,infectionK,time_step)
%end

%%%Modify above codes
for i = 1:n_tot_samples
    header = num2str(i);
    comment = ['[&' annotation '="' num2str(indiv_sampled_antigenicDrift(i),'%05.2f') '",Netcharge.set="",Bindingscr.set="' num2str(indiv_sampled_binding_s(i)) '"]'];
    nexus_names{i} = strcat('sample', int2str(i), '_' ,num2str(indiv_sampled(i)), '_', num2str(seq_times(i)));
    names{i} = [header comment];
end

% The ancestral lineages of all sampled nodes
parentLineages = GetIndividualsLineages_indiv(n_seqs, indiv_sampled, parent);

b = NaN*ones((n_seqs - 1),2);
d = NaN*ones(2*n_seqs - 1, 1);
loc_b = 1;

curr_indiv = indiv_sampled;
complete_indiv = indiv_sampled;
complete_seq_times = seq_times;

while length(curr_indiv)>0
    
    % curr_indiv: sampled individuals
    % parentLineages: parental lineages of sampled individuals   
    % coal_daughters: two external tips (vid) from sampled individuals
    % coal_parent: internal node (vid) from all individuals including
    % unsampled ones

    [coal_daughters, coal_parent, timeOfCoalescence] = FindMostRecentCoalescence_indiv(curr_indiv, parentLineages, births, deaths);
    %if coal_daughters(1) == coal_daughters(2)
    %  parent_index = find(curr_indiv==coal_parent);
    %  curr_indiv(parent_index(end)) = [];
    %  continue;
    %end
    if timeOfCoalescence == -Inf
       break;
    end
    
    indiv1_index = find(complete_indiv == coal_daughters(1));
    indiv2_index = find(complete_indiv == coal_daughters(2));

    % update distance matrix
    d(indiv1_index(end), 1) = complete_seq_times(indiv1_index(end)) - timeOfCoalescence;
    d(indiv2_index(end), 1) = complete_seq_times(indiv2_index(end)) - timeOfCoalescence;

    % update branch matrix
    if indiv1_index(end) == indiv2_index(end)
      disp('something is wrong');
    end
    b(loc_b, 1:2) = [indiv1_index(end) indiv2_index(end)];
    
    % create internal nodes
    %names{n_tot_samples + loc_b} = strcat('node', int2str(loc_b), '_', num2str(coal_parent) ,'_' ,num2str(timeOfCoalescence)); % add to name
    header = num2str(n_tot_samples+loc_b);

    %use starting binding avidity
    comment = ['[&' annotation '="' num2str(antigenicTot(coal_parent),'%05.2f') '",Netcharge.set="",Bindingscr.set="' num2str(binding(coal_parent)) '"]'];
    %nexus_names{i} = strcat('node', int2str(i), '_' ,num2str(indiv_sampled(i)), '_', num2str(seq_times(i)));
    names{n_tot_samples + loc_b} = [header comment];
    
    if coal_parent == coal_daughters(1) || coal_parent == coal_daughters(2)
      disp('one of the node is the ancestor of the other');
    end    
    
    complete_indiv = [complete_indiv coal_parent]; % add to individual list
    complete_seq_times = [complete_seq_times timeOfCoalescence]; % add to seq time list
    
    loc_b = loc_b + 1;                                                                                                                                                                                                                                                                                                                                                                                                                                          
    

    
    % erase coal_daughters from curr_indiv
    loc1_inCurr = find(curr_indiv == coal_daughters(1));
    loc2_inCurr = find(curr_indiv == coal_daughters(2));
    if coal_daughters(1)==coal_daughters(2)
      loc1_inCurr = loc1_inCurr(1);
      loc2_inCurr = loc2_inCurr(2);
    end
    if loc1_inCurr(1)>loc2_inCurr(1)
        curr_indiv(loc1_inCurr(1)) = [];
        curr_indiv(loc2_inCurr(1)) = [];
    else
        curr_indiv(loc2_inCurr(1)) = [];
        curr_indiv(loc1_inCurr(1)) = [];
    end
    
    % add coal_parent to curr_indiv list
    %%what happened if coal_parent is same as the daughters
    curr_indiv = [curr_indiv coal_parent];
    
    % redo parentLineages
    n_seqs = length(curr_indiv)
    %n_linages = length(parentLineages);
    add_parentLineages = GetIndividualsLineages_indiv(1, curr_indiv(end), parent); 
    parentLineages(end+1).lin = add_parentLineages.lin;
    parentLineages(end).tip = add_parentLineages.tip;
    if loc1_inCurr(1)>loc2_inCurr(1)
        parentLineages(loc1_inCurr) = [];
        parentLineages(loc2_inCurr) = [];
    else
        parentLineages(loc2_inCurr) = [];
        parentLineages(loc1_inCurr) = [];
    end
    %clear parentLineages
    %parentLineages = GetIndividualsLineages_indiv(n_seqs, curr_indiv, parent);        
    
end


% the remaining lineages do not coalesce
% figure out roughly the lambda
%lambda_approx = 1/min([d(indiv1_index(end), 1), d(indiv2_index(end), 1)])
%K_approx = lambda_approx/((n_seqs + 1)*n_seqs);

 
% link together arbitrarily
while 1
    if n_seqs == 1
        break;
    end
    
    perm_indiv = randperm(length(curr_indiv));
    coal_daughters = curr_indiv(perm_indiv(1:2));
    
    indiv1_index = find(complete_indiv == coal_daughters(1));
    indiv2_index = find(complete_indiv == coal_daughters(2));
    
    % erase coal_daughters from curr_indiv
    loc1_inCurr = find(curr_indiv == coal_daughters(1));
    curr_indiv(loc1_inCurr(1)) = [];
    loc2_inCurr = find(curr_indiv == coal_daughters(2));
    curr_indiv(loc2_inCurr(1)) = [];
    
    if indiv1_index(end) == indiv2_index(end)
      disp('something is wrong');
    end
    %colesce1
    %branch_length = exprnd(1/(K_approx*n_seqs*(n_seqs-1)));
    %timeOfCoalescence = min([complete_seq_times(indiv1_index) complete_seq_times(indiv2_index)]) - branch_length;
    
    %colesce2
    recovery_rate = 1/3.3; % in days
    %%%%%%num_infected = 10;  % this affect the coalescence distance
    num_infected = 10;
    lambda_val = n_seqs*(n_seqs-1)*recovery_rate/num_infected;
    branch_length = exprnd(1/lambda_val);
    timeOfCoalescence = min([complete_seq_times(indiv1_index) complete_seq_times(indiv2_index)]) - branch_length;
    
    who
    %colesce3
    %timeOfCoalescence = 0; % just have all the remaining coalesce at time t = -50
    %timeOfCoalescence = time.start-7; % just have all the remaining coalesce at time t = -50
    %timeOfCoalescence = epi_params.tRange_stoch(1)-7; 
    
    d(indiv1_index, 1) = complete_seq_times(indiv1_index) - timeOfCoalescence;
    d(indiv2_index, 1) = complete_seq_times(indiv2_index) - timeOfCoalescence;
    
    if (find(d<0))
      disp 'error: distance is negative.';
    end

    b(loc_b, 1:2) = [indiv1_index(end) indiv2_index(end)];
    % timeOfCoalescence < 0
    header = num2str(n_tot_samples+loc_b);
    comment = ['[&' annotation '="' num2str(0,'%05.2f') '",' 'Netcharge.set="",Bindingscr.set="' num2str(mean(binding(1:100))) '"]'];
    names{n_tot_samples + loc_b} = [header comment];
    loc_b = loc_b + 1;
    
    coal_parent = max(complete_indiv) + 1;
    
    complete_indiv = [complete_indiv coal_parent];
    complete_seq_times = [complete_seq_times timeOfCoalescence];
    
    
    % add coal_parent to curr_indiv list
    curr_indiv = [curr_indiv coal_parent];
    
    % redo parentLineages
    n_seqs = length(curr_indiv)
    
end
    
d(end) = 0;
% output tree array, birth, death and names
% save([folder 'tree_' num2str(n_tot_samples) '.mat'], 'b','d', 'names');

if (find(d<0))
      disp 'error: distance is negative.';
      disp 'change distance to zero if it is negative';
      d(find(d<0))=0;
end

%% resolve the tree to be fully bifurcated
sb = sort(b(:));
for (i=1:length(sb))
  if sb(i)~=i
     loss_index = find(b(:,1)==sb(i)); 
     u = unique(sb);
     un = histc(sb,u);
     rep = sb(un>1);
     dup_index = find(b(:,1)==rep(1));
     b(dup_index,1) = i;
     rep(1) = [];
     if isempty(rep)
         break;
     else
         sb = sort(b(:));
         continue;
     end
  end
end

tree = phytree(b,d, names);
%view(tree)
plot(tree);
pos_s = regexp(filename,'\.');
save([filename '_simple.mat'],'b','d','names');
saveas(gcf,[filename_infectionTreeData '.jpg']);
if epi_params.display==0
  close(gcf);
end

phytreewrite(outfile_tree, tree, 'BranchNames', false);

if epi_params.savetree == 1 | epi_params.savetree == 2

% write nexus tree
pos_s = regexp(outfile_tree,'\.');
nexusoutfile = [outfile_tree(1:pos_s-1) '.nex'];

fid=fopen(nexusoutfile, 'w');
seq = ['#NEXUS\n\n' 'Begin taxa;\n' 'Dimensions ntax=' num2str(n_tot_samples) ';\n' 'Taxlabels\n'];
fprintf(fid, seq);
for i=1:length(nexus_names)
  fprintf(fid, [nexus_names{i} '\n']);
end
seq = [';\n' 'End;\n\n' 'Begin trees;\n' 'Translate\n'];
fprintf(fid, seq);
for i=1:length(nexus_names)
    if i<length(nexus_names)
      fprintf(fid, [num2str(i) ' ' nexus_names{i} ',\n']);
    else
      fprintf(fid, [num2str(i) ' ' nexus_names{i} '\n']);
    end
end
seq = [';\n' 'tree TREE1 = [&R]'];
fprintf(fid, seq);

% append phytree into nexus file 
%tree_str=fileread(treeoutfile);
tree_str = fileread(outfile_tree);
expression = '\';
replace = '';
new_tree_str = regexprep(tree_str,expression,replace);
fprintf(fid, new_tree_str);
fprintf(fid, ['\nEnd;\n']);
fclose all;
end

if epi_params.savetree == 2
    
    % output tree data
%save(outfile_treeData); %save time, don't save output now
tree_str = fileread(outfile_tree);
expression = '\';
replace = '';
new_tree_str = regexprep(tree_str,expression,replace);
fid=fopen(outfile_tree, 'w');
fprintf(fid, new_tree_str);
fclose all;

% output internal and external nodes lengths
pos_s = regexp(outfile_tree,'\.');
outfile_tree_allnodes = [outfile_tree(1:pos_s-1) '.allnodes' outfile_tree(pos_s:end)];
% write phytree
phytreewrite(outfile_tree_allnodes, tree, 'BranchNames', true);
tree_str = fileread(outfile_tree_allnodes);
expression = '\';
replace = '';
new_tree_str = regexprep(tree_str,expression,replace);
fid=fopen(outfile_tree_allnodes, 'w');
fprintf(fid, new_tree_str);
fclose all;

% write nexus tree for all nodes
pos_s = regexp(outfile_tree,'\.');
nexusoutfile = [outfile_tree(1:pos_s-1) '_allnodes' '.nex'];
fid=fopen(nexusoutfile, 'w');
seq = ['#NEXUS\n\n' 'Begin taxa;\n' 'Dimensions ntax=' num2str(n_tot_samples) ';\n' 'Taxlabels\n'];
fprintf(fid, seq);
for i=1:length(nexus_names)
  fprintf(fid, [nexus_names{i} '\n']);
end
seq = [';\n' 'End;\n\n' 'Begin trees;\n' 'Translate\n'];
fprintf(fid, seq);
for i=1:length(nexus_names)
    if i<length(nexus_names)
      fprintf(fid, [num2str(i) ' ' nexus_names{i} ',\n']);
    else
      fprintf(fid, [num2str(i) ' ' nexus_names{i} '\n']);
    end
end
seq = [';\n' 'tree TREE1 = [&R]'];
fprintf(fid, seq);
% append phytree into nexus file 
%tree_str=fileread(treeoutfile);
tree_str = fileread(outfile_tree_allnodes);
fprintf(fid, tree_str);
fprintf(fid, ['\nEnd;\n']);
fclose all;
%delete(outfile_tree_allnodes); % keep the nexus but delete the original allnodes tree file 
%disp(['delete the temporary file: ' outfile_tree_allnodes]);
end

delete(outfile_tree);
disp(['delete the temporary file: ' outfile_tree]);
end

function v = getBindingAvidity(t,v,k)
    p = path
    path(p,'lib/');
    if ~exist('params')
        params.p = 4;
        params.r = 70;
        params.a = 0.7;
        params.b = 3;
        params.c = 0.5;
    end
    t0 = 0;
    dT = t;
    nsteps = 30;
    tspan = linspace(t0,dT,nsteps);
    pars.k = k;
    [v1] = ode2(@odef_v_change, tspan, v, pars);
    v = v1(end);
end



