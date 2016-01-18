%Purpose: calculate the number of each amino acid occured in given sequences
%Parameters: None
%Input files: FASTA files downloaded from Influenza Viruses Database
%Output: 
%   1) aa_array: cell array of #amino acids for each fasta sequence
%   2) aa_charge: isolation year and net charges for each fasta sequence
%   3) aa_charge_sum: year, net charges 
function [aa_array aa_charge aa_charge_sum] = group_aa_dist_byngs_byK(viruses)

%1)iso_year, 2)charge, 3)ngs and 4)number
aa_array = [];
ngs_index = 4;
iso_index = 2; %isolation date
charge_index = 3; %netcharge
K_index = 6; % #of previous infection

ac = []; %year, charge, ngs, number
%for i=1:length(viruses(:,1))
    v = [viruses(:,K_index) viruses(:,ngs_index) viruses(:,charge_index)];
    ac = push_data(v,ac);
%end
aa_array = ac;
disp 'done';
function plot_aa_charge(ac)
  siz = length(ac);
  %figure;
  %hold on;
  for i=1:siz
      yr = round(ac(i,1));
      charge = ac(i,2);
      amount = ac(i,3);
      if(amount<5)
        plot(yr,charge,'g.', 'MarkerSize', 5);
      elseif(amount<7)
        plot(yr,charge,'g.', 'MarkerSize', 10);
      elseif(amount<15)
           plot(yr,charge,'g.', 'MarkerSize', 15);
      else
          plot(yr,charge,'g.', 'MarkerSize', 20);
      end
  hold on;
  end
end


function [ac]= push_data(v,ac)

    for i=1:length(v(:,1))
        infK=round(v(i,1));
        ngs=v(i,2);
        charge=v(i,3);
        %d(i,1:5), 1=year, 2=#ngs, 3=netcharge
        %%4=#pos aas, 5=#neg aas, 6=#pos+#neg
        if isempty(ac)
         ac = [infK ngs charge 1];
         continue; 
        end
        %loc = find(ac(:,1)==yr & ac(:,3)==charge);
        loc = find(ac(:,1)==infK & ac(:,3)==charge &ac(:,2)==ngs);
        
        if(isempty(loc))
            ac = [ac; infK ngs charge 1];
        else
            ac(loc,4) = ac(loc,4) + 1;
        end
    end
    %ac(1,:) = [];
end

end


