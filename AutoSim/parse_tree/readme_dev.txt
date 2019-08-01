Aug 02, 2013
proj = 'ancestor_20130713';
annotate_tree_binding(['dat/' proj '/hm_h3n2_ny_dna_beast_1993.nw.mcc.(time).trees'],['dat/' proj '/hm_h3n2_ny_charge_1993_2006.mat'],['dat/' proj '/hm_h3n2_ny_binding_1993.mat'])

proj = 'test_burnin-20130731';
[tree, n, b, nm, taxa] = read_newick_seq(['dat/' proj '/hm_h1n1_noram_1995_2006_nogap_bur95.nw.mcc.(time).trees']);
save(['dat/' proj '/newick_elements.mat'], 'b', 'n', 'nm', 'tree', 'taxa');
taxa2fasta(taxa); ->ancestral_sequences.fas

Jul 31, 2013 (2)
Substitute insertion --- by CAA.
Rerun BEAST in folder nogap-20130731
proj = 'nogap-20130731';
[tree, n, b, nm, taxa] = read_newick_seq(['dat/' proj '/hm_h1n1_noram_1995_2006_nogap.nw.mcc.(time).trees']);
save(['dat/' proj '/newick_elements.mat'], 'b', 'n', 'nm', 'tree', 'taxa');
taxa2fasta_new(taxa); ->ancestral_sequences.fas  % consensus sequence of posterior distribution of sequences
Use MEGA to transfer ancestral_sequences.fas -> ancestral_sequences.prot.fas
main_aa_dist([proj '/ancestral_sequences.prot.fas'],'hm_h1n1_noram_charge_1995_2006_nogap.mat','H1N1')
change directory to noram9506-20130730
annotate_tree('dat/noram9506-20130730/hm_h1n1_noram_1995_2006.nw.mcc.(time).trees',['dat/' proj '/hm_h1n1_noram_charge_1995_2006_nogap.mat'])
Insert nexus annotation back, e.g, 'tree TREE1 = [&R]'.

Jul 31, 2013 (1)
Annotate viral netcharge on H1N1 phylogenetic tree
proj = 'noram9506-20130730';
[tree, n, b, nm, taxa] = read_newick_seq('dat/noram9506-20130730/hm_h1n1_noram_1995_2006.nw.mcc.(time).trees');
save('dat/noram9506-20130730/newick_elements.mat', 'b', 'n', 'nm', 'tree', 'taxa');
taxa2fasta(taxa);
Use MEGA to transfer ancestral_sequences.fas -> ancestral_sequences.prot.fas
main_aa_dist([proj '/ancestral_sequences.prot.fas'],'hm_h1n1_noram_charge_1995_2006.mat','H1N1') -> hm_h3n2_ny_charge_1993 (typo, should be h1n1)
annotate_tree('dat/noram9506-20130730/hm_h1n1_noram_1995_2006.nw.mcc.(time).trees','dat/noram9506-20130730/hm_h1n1_noram_charge_1995_2006.mat')
Insert nexus annotation back, e.g, 'tree TREE1 = [&R]'.


Jul 2, 2013
Target template: hm_h3n2_mcc_trait3.nx.trees
Can show different colors for different traits on branches.

Flows:
1.  Open nexus files and save only newick data into a separate newick file.
2.	Read the newick file with sequences. (read_newick_seq.m -> save tree topology and taxa)
3.	Transfer sequences in taxa to fasta. (taxa2fasta.m -> ancestral_sequences.fas)
4.	Convert into translated sequences using MEGA. (ancestral_sequences.fas-> ancestral_sequences_prot.fas)
5.	Calculate netcharge. (main_aa_dist(infile, outfile, 'H3N2') -> hm_h3n2_ny_charge_1993)
6.	Annotate newick tree with netcharge. (annotate_tree.m -> hm_h3n2_ny_dna_beast_1993.traits.nw.mcc.(time).trees)
6.1 Annotate newick tree with netcharge and binding score. (annotate_tree_binding.m)
7.	Copy the nexus header into the newick file. (hm_h3n2_mcc.traits.nx.trees)


proj = 'ancestor_20130702'
Steps for annotate viral netcharge on phylogenetic tree
>> [tree, n, b, nm, taxa] = read_newick_seq('dat/ancestor_20130702/hm_h3n2_ny_dna_beast_1993_2006.nw.mcc.(time).trees')
>> save('newick_elements.mat', 'b', 'n', 'nm', 'tree', 'taxa')
>> taxa2fasta(taxa)
Use MEGA to transfer ancestral_sequences.fas -> ancestral_sequences.prot.fas
>> main_aa_dist([proj '/ancestral_sequences.prot.fas'],'hm_h3n2_ny_charge_1993_2006.mat','H3N2') -> hm_h3n2_ny_charge_1993
>> annotate_tree(['dat/' proj '/hm_h3n2_ny_dna_beast_1993.nw.mcc.(time).trees'],['dat/' proj '/hm_h3n2_ny_charge_1993_2006.mat'])
>> annotate_tree_binding('dat/ancestor_20130713/hm_h3n2_ny_dna_beast_1993.nw.mcc.(time).trees','dat/ancestor_20130713/hm_h3n2_ny_charge_1993.mat','dat/ancestor_20130713/hm_h3n2_ny_binding_1993.mat')


%Newick to pairs
[tree, n, b, nm, taxa] = read_newick_seq('dat/ancestor_20130702/hm_h3n2_ny_dna_beast_1993_2006.nw.mcc.(time).trees')
-> save tree topology as hm_h3n2_ny_dna_beast_1993_2006.topo.nw.mcc.(time).trees
pairs = string2pairs('dat/ancestor_20130702/hm_h3n2_ny_dna_beast_1993_2006.topo.nw.mcc.(time).trees',665);
%618 ABJ53482_CY016995_2006/04/06_11_17,
[anc]= getParentList(pairs, 618); %identify the main trunk



%Pairs to Newick
tree = pairs2string(pairs,b); %need branch lengths
fileID = fopen(outfile,'w');
fprintf(fileID,'%s',tr);
fclose(fileID);

%Annotate traits on Newick
annotate_tree('tree2.txt','dat/ancestor_20130702/hm_h3n2_ny_charge_1993_2006.mat');