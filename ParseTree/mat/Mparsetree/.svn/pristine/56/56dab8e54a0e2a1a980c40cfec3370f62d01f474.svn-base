function tree = pairs2string(pairs, b)
% function pairs = string2pairs(tree)
% Converts string record of tree to a pair one.
%
% PHYLLAB toolbox v1.1
%need to read nm, otherwise, the id is not right

[k,l] = size(pairs); %k=#row, l=#col
n = k/2; 
%tr = tree;
tr = '';

posO = 1;
flagO=0;

% parsing the input
tmp_pairs = pairs; 
global node_array;
global branch;
branch = b;
node_array = [];
m = max(tmp_pairs);
while  m>n 
    %editing by HY
    %ch -> children
    if(isnumeric(m))
    else
        m = str2num(m);
    end
    ch = find(tmp_pairs == m);
    %if only 1 child, then S->C
    if (length(ch)==1)
        tr = sprintf('%s',['(' num2str(ch) ':0)']);
        m = ch;
    end
    %if 2 children, then S-> (C1,C2)
    if (length(ch)==2)
        tree = parse(tr,m, tmp_pairs)
        return;
    end;
end;

end


function tree = parse(tree,m, tmp_pairs)
    global node_array; 
    global branch;
    ch = find(tmp_pairs == m);
    b_length = branch(ch);

    if (length(ch)==1)
        
        if(ch(1)<10)
                ch_str = ['00' int2str(ch(1))];
        else
                ch_str = ['0' int2str(ch(1))];
        end
        
        %mloc: the location of the symbol you want to replace
        %mlength: the liegth of the words 
        m = [num2str(m) ':'];
        mloc = find_loc(tree,num2str(m))
        mlength = length(num2str(m)) - 1; 
        %add branch length
        %if (mloc ~= -1)
        %    tree = sprintf('%s%s%d%s%d%s%s', [tree(1:mloc-1) '(' num2str(ch(1)) ',' num2str(ch(2)) ')' tree(mloc+1:length(tree))]);
        %end
        if (mloc ~= -1)
            %replace the symbol
            tree = sprintf('%s', [tree(1:mloc-1) '(' ch_str(1) ':' num2str(b_length) ')' tree(mloc+mlength:length(tree))]);
        end
        %push children nodes into array
        push_array(ch(1));
    end;
    %if 2 children, then S-> (C1,C2)
    if (length(ch)==2)
        ch_str1 = '';
        ch_str2 = '';
        if(ch(1)<10)
            ch_str1 = ['0' int2str(ch(1))];
        else
            ch_str1 = num2str(ch(1));
        end
        if(ch(2)<10)
            ch_str2 = ['0' int2str(ch(2))];
        else 
            ch_str2 = num2str(ch(2));
        end
        %mloc is the location of the symbol you want to replace 
        %test: avoid search duplicate result
        m = [num2str(m) ':'];
        mloc = find_loc(tree,num2str(m))
        mlength = length(num2str(m)) - 1; 
        %add branch length
        %if (mloc ~= -1)
        %    tree = sprintf('%s%s%d%s%d%s%s', [tree(1:mloc-1) '(' num2str(ch(1)) ',' num2str(ch(2)) ')' tree(mloc+1:length(tree))]);
        %end
        if (mloc ~= -1)
            %replace the symbol
            %tree = sprintf('%s%s%d%s%d%s%s', [tree(1:mloc-1) '(' num2str(ch(1)) ',' num2str(ch(2)) ')' tree(mloc+mlength:length(tree))]);
            tree = sprintf('%s', [tree(1:mloc-1) '(' ch_str1 ':' num2str(b_length(1)) ',' ch_str2 ':' num2str(b_length(2)) ')' tree(mloc+mlength:length(tree))]);
        end
        %push children nodes into array
        push_array(ch(1));
        push_array(ch(2));
    end;
    %then find the nodes to process repeatedly
    tmp_ch = pop_array();
    if(tmp_ch == 52)
        disp(tmp_ch);
    end
    if(isempty(tmp_ch))
      return;
    end
    tree = parse(tree, tmp_ch, tmp_pairs);
end


function loc = find_loc(str, pattern)            % Subfunction
  %pat = [pattern ',']; 
  
  loc = strfind(str, pattern)
  %disp(loc)
  if(isempty(loc))
      loc = -1;
  end
end  
% Push every nodes into array when children found
function push_array(c)
  global node_array;
  siz1 = length(node_array);
  node_array(siz1+1) = c;
end

function a = pop_array()
  global node_array;
  siz2 = length(node_array);
  if(siz2 == 0)
      a = [];
      return;
  end
  a = node_array(siz2);
  node_array(siz2) = [];
end