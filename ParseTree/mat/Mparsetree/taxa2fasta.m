function [] = taxa2fasta( taxa, outdir)
%function [] = taxa2fasta(taxa)
% Detailed explanation goes here
% Save all the external sequences and internal reconstructed sequences into
% FASTA file.
% output: ancestral_sequences.fas
% taxa struct:
% taxa.id = sequential ID in newick file
% taxa.name = original ID in nexus tree file
% taxa.annotation = annotation in [...]

if ~exist('outdir')
  outdir = './';
end

for i =1:length(taxa)
   if i==189
     disp check;
   end
   Header(i) = {num2str(taxa(i).id)};
   str = taxa(i).annotation;
   pos_s = regexp(str, '="[AaGgTtCc]')+2; %start at "AGTC...
   %pos_e = regexp(str, '[AaGgTtCc]",'); %stop at ",
   %pos_e = regexp(str, '[AaGgTtCc^"]"(,|\])'); %stop at ", or "]  Beware A","A in set
   %pos_l = regexp(str(pos_s:end), '[AaGgTtCc^"]"(,|\])');
   pos_l = regexp(str(pos_s:end), '[AaGgTtCc^"](+|",|"\])'); %stop at ", or "] or + Beware A","A in set
   if (length(pos_l)>1)
     %disp 'warning: multiple hits.';
     pos_l = pos_l(1);
   end
   pos_e = pos_s+pos_l-1;  
   
   %pos_s = regexp(str, 'set={"[AaGgTtCc]')+6;
   %pos_e = regexp(str, '[AaGgTtCc]"}');
   leng_std = pos_e-pos_s;
   if pos_e-pos_s ~= leng_std
     disp 'error: sequences length not consistent';
   end
   
   if (pos_s(1)>0 & pos_e(1)>pos_s(1))
       Sequence(i) = {char(str(pos_s:pos_e))};
   else
       Sequence(i) = '';
   end
end
   if (strcmp(outdir(end), '/'))
      fastawrite([outdir 'ancestral_sequences.fas'], Header, Sequence);
   else 
      fastawrite([outdir '/ancestral_sequences.fas'], Header, Sequence);
   end
   
   disp gofasta;
end

