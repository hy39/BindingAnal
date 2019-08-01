function pairs = string2pairs(file,n)
% function pairs = string2pairs(tree)
% Converts string record of tree to a pair one.
%
% PHYLLAB toolbox v1.1


ff = fopen(file,'r');
if ff == -1 
   n=-1;	tree = [];	
   disp('File open error in read_newick.');
   return;
end;
% The tree can be recorded in multiple lines. Removing LF, CR, etc...
[tree, k] = fscanf(ff,'%c'); 


[k,l] = size(tree);
tr = tree;
pairs = zeros(n*2-1,1); % format from the VOSTORG package; add zero to the last element

anc = n+1;
posO = 1;
flagO=0;

% parsing the input
i=1;
while  i<=l 
    if tr(i) == '(' 
        flagO = 1;posO = i;
    end;
    if tr(i) ~= ')' & tr(i) ~= '(' & tr(i) ~= ','
        flagO = 2;
    end
    if abs(tr(i))<abs(' ')
        break; 
    end;
    if tr(i) == ')' & flagO == 2
		% now we process the nest
      if anc==2*n break; end;
		sbst = tr(posO:i);
		[pr, np] = sscanf(sbst,'(%d,%d)');
        %[pr, np] = sscanf(sbst,'(%s,%s)'); %test if I scan string
      pairs(pr(1)) = anc;         
      pairs(pr(2)) = anc;                       
      % now we rewrite the tree string
      tr = sprintf('%s%d',tr(1:posO-1),anc,tr(i+1:l));
      anc=anc+1;
      i=0;
      [k,l]=size(tr);
   end;
	i=i+1;     
end;
%pairs(anc-1)=anc;   
