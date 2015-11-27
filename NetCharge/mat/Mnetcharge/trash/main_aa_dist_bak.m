%calculate the number of each amino acid occured in given sequences

function [aa_array aa_dist aa_charge] = main_aa_dist()

file_name = 'protein_ha.fa';
[Header, Sequence] = fastaread(['dat/' file_name]);
siz = length(Header);
%aa_array;
aa_charge = [];

for i=1:siz
    aa = get_aa_num(char(Sequence(i)));
    aa_array(i,1) = aa;
    aa_charge(i,1) = str2num(get_year(char(Header(i))));
    aa_charge(i,2) = aa.R + aa.H + aa.K - aa.D - aa.E;
    if(i == 1)
            aa_dist.A = aa.A;
            aa_dist.R = aa.R;
            aa_dist.N = aa.N;
            aa_dist.D = aa.D;
            aa_dist.C = aa.C;
    end
    
    aa_dist.A = aa_dist.A + aa.A;
    aa_dist.R = aa_dist.R + aa.R;
    aa_dist.N = aa_dist.N + aa.N;
    aa_dist.D = aa_dist.D + aa.D;
    aa_dist.C = aa_dist.C + aa.C;
%    aa_dist.A = aa.Q;
%    aa_dist.A = aa.E;
%    aa_dist.A = aa.G;
%    aa_dist.A = aa.H;
%    aa_dist.A = aa.I;
%    aa_dist.A = aa.L;
%    aa_dist.A = aa.K;
%    aa_dist.A = aa.M;
%    aa_dist.A = aa.F;
%    aa_dist.A = aa.P;
%    aa_dist.A = aa.S;
%    aa_dist.A = aa.T; 
%    aa_dist.A = aa.W;
%    aa_dist.A = aa.Y;
%    aa_dist.A = aa.V;
%    aa_dist.A = aa.B;
%    aa_dist.A = aa.Z; 
%    aa_dist.A = aa.X; 

end




    function yr = get_year(h)
        [acc r]= strtok(h, '_');
        [yr r]= strtok(r, '_');
    end
end


