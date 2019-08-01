%Subplot1 f(k,V)
%Subplot2 g(V) = exp(-a*V.^b)
%Subplot3 Beta(k,V)
%Subplot4 dBeta(k,V)/dV
%Author: Hsiang-Yu Yuan
%1st version: Jul 16, 2012
%Rho = f x g
%2nd verison: Aug 01, 2013
%Rho = 1 - R0^-1
%R0 = f x g x n
function [ Vp ] = getVp(MaxK )
%Binding avidity range
V = 0:0.01:2;

%Transmission parameters
p = 4;
%r = 1;
r = 70;
b = 3;
a = 0.7;
c = 0.5; % contact rate
nv = 4; % average copies number of each virion

%Epidemiological parameters
mu_val = 1/(70*365.25); %lifespan = 70y
gamma_val=1/5; %infecious period should change to 5d


Trans_array = []; %Array of probability of evasion by immune system
Trans_Pr_array = [];

N_Reinfect = MaxK+1; 

figure;
%colorm = {[1 1 0];[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1];[0 0 0]}; 
%

for k = 0:N_Reinfect-1
%for k=1:20:21;   
   P_Ab = exp(-p*(V+1));
   
   %%%P_Trans: f(k,V)
   %%%P_Trans_Pr: f'(k,V)
   P_Trans = (1-P_Ab).^(r*k); 
   if k>=1
       P_Trans_Pr = r*k*p.*((1-P_Ab).^(r*k-1)).*(P_Ab);
   elseif k == 0
       P_Trans_Pr = zeros(1,length(V));
   end
   %subplot(2,2,1); plot(V,P_Trans,'color',cell2mat(colorm(k+1))); hold on;
   Trans_array(k+1,:) = P_Trans;
   Trans_Pr_array(k+1,:) = P_Trans_Pr;
end
P_Rep = exp(-a*V.^b);
P_Rep_Pr = -a*b*exp(-a*V.^b).*(V.^(b-1));

for i=1:length(Trans_Pr_array(:,1))
%R0 = f(k,V)g(V)n
%rho = 1 - 1/R0
%beta = c x rho
R0_Trans = Trans_array(i,:).*P_Rep.*nv;
%find Vp that produce the maximum R0 
v_elem = find(R0_Trans==max(R0_Trans));
Vp(i) = V(v_elem);
Rho_Trans = 1 - R0_Trans.^-1;
Rho_Trans(find(Rho_Trans<0))=0;
B_Trans = c.*Rho_Trans;
end
   

