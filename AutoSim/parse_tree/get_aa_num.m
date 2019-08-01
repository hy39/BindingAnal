function [L] = get_aa_num(seq)

  L.A = length(find(seq=='A'));%      Alanine
  L.R = length(find(seq=='R'));%      Arginine
  L.N = length(find(seq=='N'));%      Asparagine
  L.D = length(find(seq=='D'));%      Aspartic acid (Aspartate)
  L.C = length(find(seq=='C'))%       Cysteine
  L.Q = length(find(seq=='Q'))
  L.E = length(find(seq=='E'))
  L.G = length(find(seq=='G'))
  L.H = length(find(seq=='H'))
  L.I = length(find(seq=='I'))
  L.L = length(find(seq=='L'))
  L.K = length(find(seq=='K'))
  L.M = length(find(seq=='M'))
  L.F = length(find(seq=='F'))
  L.P = length(find(seq=='P'))
  L.S = length(find(seq=='S'))
  L.T = length(find(seq=='T'))
  L.W = length(find(seq=='W'))
  L.Y = length(find(seq=='Y'))
  L.V = length(find(seq=='V'))
  L.B = length(find(seq=='B'))
  L.Z = length(find(seq=='Z'))
  L.X = length(find(seq=='X'))
  
  L.NGS = length(regexp(seq, 'N[^P][S|T]'));
  L.locs = regexp(seq, 'N[^P][S|T]')';
  
%                A       Alanine
%                R+      Arginine
%                N       Asparagine
%        Asp     D-      Aspartic acid (Aspartate)
%        Cys     C       Cysteine
%        Gln     Q       Glutamine
%        Glu     E-      Glutamic acid (Glutamate)
%        Gly     G       Glycine



%        His     H+      Histidine
%        Ile     I       Isoleucine
%        Leu     L       Leucine
%        Lys     K+      Lysine
%        Met     M       Methionine




%        Phe     F       Phenylalanine
%        Pro     P       Proline
%        Ser     S       Serine
%        Thr     T       Threonine
%        Trp     W       Tryptophan



%        Tyr     Y       Tyrosine
%        Val     V       Valine
%        Asx     B       Aspartic acid or Asparagine
%        Glx     Z       Glutamine or Glutamic acid.
%        Xaa     X       Any amino acid.
end
