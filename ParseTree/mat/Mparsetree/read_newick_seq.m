function [t, n, b, nm, taxa] = read_newick_seq(file)
% function [tree,total,branches,names] = read_newick(tstr)
% Reads tree from standart file in Newick format 
% 
% Description of this format can be found at
% http://evolution.genetics.washington.edu/phylip/newicktree.html
%
% n -- the number of OTUs
% tree -- a string with the tree
% branches -- branch length
% names -- OTU names
% If no branch length provided, thay are set to ones.
%
%
% PHYLLAB toolbox v1.1
% ex. [tree n b nm, taxa] = read_newick_seq('hm_h3n2_mcc.nx.trees')

% opening file
ff = fopen(file,'r');
if ff == -1 
   n=-1;	tree = [];	
   disp('File open error in read_newick.');
   return;
end;

% The tree can be recorded in multiple lines. Removing LF, CR, etc...
[tr, k] = fscanf(ff,'%c'); 
n=0;
i=1;
while i<k 
    if tr(i) <= ' '
       tr(i) =[]; 
       k=k-1;
    else
       i=i+1;
    end
    if tr(i) == ')'  
        n=n+1;
    end
end
n=n+1; %  TEMP!!!!!!!!!!!!!!!!!!!!!!!!!

% Now replacing all symbolic names to numbers!
% Even if the name is already number it will be replaced
posO = 1; % position of last '(' - to determin the start of the nest 
flagO=0; % the status of previous symbol: 1 - if punctuation, 2 - other
flagA=0; % the status of previous symbol: 1 - if annotation, 0 - normal
nm = char(n,n);
taxa = struct;
cur_tip=1;
cur_inter_node=n+1;

[l,k] = size(tr);
i=1;
tr

pos_annotation_s = 0;
pos_annotation_e = 0;
annotation = '';
while i<k
    if i>length(tr)
      break;
    end
    %skip the processes if tr(i) is in the bracket [...] Jul 1, 2013
    if (tr(i) == '[' & flagA ==0)
       pos_annotation_s = i;
       flagA = 1;
       i=i+1;
       continue;
    end

    if (tr(i) == ']' & flagA ==1)
       pos_annotation_e = i;
       flagA = 0;
       i=i+1;
       continue;
    end
    %if (pos_annotation_s > 1 && flagA == 0)
      %annotation = tr(pos_annotation_s:pos_annotation_e);
    %end
    if (flagA == 0)
    %parse internal nodes
    if (tr(i) == ':' | tr(i) == ';' ) & flagO == 1
        %cur_inter_node %disp 
        taxa(cur_inter_node).id = cur_inter_node;
        taxa(cur_inter_node).name = '';
        taxa(cur_inter_node).annotation = tr(pos_annotation_s:pos_annotation_e);
        %tr(posO+1:i-1);
        
        %if cur_inter_node > 1 nm=strvcat(nm,tr(posO+1:i-1));
        %  else nm=tr(posO+1:i-1);
        %end;
        %Add internal nodes
        %tr = sprintf('%s%d%s',tr(1:posO),cur_inter_node,tr(i:end)); 
        %[j,i]=size(sprintf('%s%d',tr(1:posO),cur_inter_node));
        %i=i+1;
        
        tr = sprintf('%s%s',tr(1:posO),tr(i:end));
        [j,i]=size(sprintf('%s',tr(1:posO)));
        i=i+1;
        cur_inter_node=cur_inter_node+1;
    end
    
	%if (tr(i) == ')' | tr(i) == '(' | tr(i) == ',' | tr(i) == ':' ) & flagO == 2
    if (tr(i) == ')' | tr(i) == '(' | tr(i) == ',' | tr(i) == ':') & flagO == 2
      cur_tip %disp 
      taxa(cur_tip).id = cur_tip;
      if pos_annotation_s~=0
      taxa(cur_tip).name = tr(posO+1:pos_annotation_s-1);
      taxa(cur_tip).annotation = tr(pos_annotation_s:pos_annotation_e);
      else
      taxa(cur_tip).name = tr(posO+1:i-1);
      end
      if cur_tip > 1 
          nm=strvcat(nm,tr(posO+1:i-1));
          else nm=tr(posO+1:i-1);
      end;
      % Assign sequential ID
      %tr = sprintf('%s%d%s',tr(1:posO),cur_tip,tr(i:end));
      %[j,i]=size(sprintf('%s%d',tr(1:posO),cur_tip));
      % Use original ID
      tr = sprintf('%s%s%s',tr(1:posO),taxa(cur_tip).name,tr(i:end));
      [j,i]=size(sprintf('%s%s',tr(1:posO),taxa(cur_tip).name));
      i=i+1;
      [j,k]=size(tr);
      cur_tip=cur_tip+1;
      flagO = 3;%%% added by HY
   end;
   if (tr(i) == '(' | tr(i) == ')' | tr(i) == ',')
        flagO = 1;posO = i;
   end;
   if tr(i) == ':' 
      flagO = 3;
   end;
   if (tr(i) ~= ')' & tr(i) ~= '(' & tr(i) ~= ',' & tr(i) ~= ';' & tr(i) ~= ':' & flagO ~= 3)
      flagO = 2;
   end;
   end
	i=i+1;  
end;

% Save trees file with viral trait
outfile = 'test.tree';
fileID = fopen(outfile,'w');
fprintf(fileID,'%s',tr);
fclose(fileID);

%***************
  % disp(tr);
%***************

% now we will deal with branch length 

[k,l] = size(tr);%return size of matrix k x l
tree=tr; % save the tree
%Create array of all zeros
%b=zeros(n,1);
b=zeros(2*n-1,1); %TEMP!!!!!!!!!!!!!!!!
anc = n+1;
posO = 1;
flagO=0;

% parsing the input
i=1;
while  i<=l 
    %tr stores newwick string
    if tr(i) == '(' 
        flagO = 1;posO = i;
    end;
    %tr(i) not equals those symbol
    if tr(i) ~= ')' & tr(i) ~= '(' & tr(i) ~= ','
       flagO = 2;
    end;
    
    if abs(tr(i))<abs(' ')
        break; 
    end;
    
    if tr(i) == ')' & flagO == 2
		% now we process the nest
      if anc==2*n break; end;
      %disp(tr);%hy
      sbst = tr(posO:i);
      
      % split string by the ',' symbol      
      [token,rem] = strtok(sbst(2:end-1),',');
		rem=rem(2:end);      
      
      % is there a ':' symbol before the ','?
      [t,r] = strtok(token,':');r=r(2:end);
      pr=str2num(t);

      if isempty(r)
         b(pr)=1;
      else
%hy
      format long;
%                      disp(r);
%                      disp(str2num(r));
         b(pr)=str2num(r);
         %b(pr)=str2num(r)*1000; %August 11, 2008 HY
      end;
      
 %     disp('Here');
      [t,r] = strtok(rem,':');r=r(2:end);
   %   disp(rem);
   %   disp(t);
   %   disp(r);
      vt=str2num(t);
      
      %if(vt==504)
      %    disp(vt);
      %end;
      
      if isempty(r)
         b(vt)=1;
      else         
         b(vt)=str2num(r);
         %b(vt)=str2num(r)*1000; %August 11, 2008 HY
      end;
      
      % now we rewrite the tree string
      % continue processing new tree string
      tr = sprintf('%s%d',tr(1:posO-1),anc,tr(i+1:l));
      %tr = sprintf('%s%s',tr(1:posO-1),tr(i+1:l));
      %disp(tr);
      %*********************
        %  disp(tr);
      %********************
      
      anc=anc+1;
      i=0;
      [k,l]=size(tr);
   end;
	i=i+1;     
end;%end of while

% remove branch length from the tree

[k,l] = size(tree);
f = 1;
j=1;  
for i=1:l
   if tree(i) == ':' f = 0; end;
   if tree(i) == ')' | tree(i) == '(' | tree(i) == ','  f = 1; end;
   if f t(j)=tree(i); j=j+1; end;
end;

% Save tree file
outfile = [file(1:regexp(file, '\.')-1) '.topo' file(regexp(file, '\.'):end)];
fileID = fopen(outfile,'w');
fprintf(fileID,'%s',t);
fclose(fileID);

