function [b] = sub_parse_branch_length(file)
% function [b] = parse_brach_length(file)
% Read the tree from Newick format 
% 
% Description of this format can be found at
% http://evolution.genetics.washington.edu/phylip/newicktree.html
%
% b -- branch length
%
% PHYLLAB toolbox v1.1
% ex. b = read_newick_seq('trees1.txt')


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
else
    


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
