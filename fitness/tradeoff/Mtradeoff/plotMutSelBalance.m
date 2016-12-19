% plot the selection coefficient with mutational effects under different
% Binding Scenario
% Fixed binding
% 

function [  ] = plotMutSelBalance()

infile = 'dat/lowV0/hostKs_low_adaptive.csv';
M = csvread( infile ) 
%V = csvread('dat/lowV0/voutput1_low_adaptive.csv',2);

fixedV = 0.8;
%MaxK = 11;

%vk = getVp(MaxK)

%X1 = normpdf(1:MaxK+1,1);
%X2 = normpdf(1:MaxK+1,5);
%X3 = normpdf(1:MaxK+1,8);
delta = 1:15;

%X = X1;
%X0 = 1-sum(X);
%X(1) = X(1) + X0;
%sk = X;
v1 = 0;
v2 = 1;
%del = 0; %no antigenic change
del = 2; %antigenic change 2 as our simulation
%f = @(v)-getOptV( v, del, sk );
%vopt = fminunc(f,v1,v2);

figure;
hold on;
ti=0;
tii = 0;
totalN = sum(M(1,:));
T = [1 200];;
for tid=1:length(T) 
    %-- draw R curve only every 5 time units
    %test 1 200 
    t = T(tid);
    Sk = M(t,1:50);
    
    X = Sk;
    %X = X./sum(X);
    X = X./totalN;
    sk = X;
    f = @(v)-getOptV( v, del, sk );
    Vopt = fminunc(f,v1,v2);
    Vnaive = 0;
    Ropt = getPopR( del, sk, Vopt);
    Rnaive = getPopR( del, sk, Vnaive);

    
    

    
    %%----------------------------
   
    V = 0:0.05:1.5;
    for i=1:length(V)
        v = V(i);

        %R = getPopR( del, sk, vopt);
        R(i) = getPopR( del, sk, v);
        %Rin(i) = getPopRin( del, sk, v);
        Rpop(i) = getPopR( del, sk, v);
        Inewt = getBeta( sk, del, v ).*sk; % estimate the newly infected individuals
        Rin(i) = getPopRin( del, Inewt, v);
    end
    
    %plot3(V, -t*5*(repmat(1,length(V))), R, 'Color',[0.4,0.4,0.4]);
    plot(V, Rin/max(Rin));
    hold on;
    plot(V, R/max(R));
    Vfinal = 0.47;
    Vini = 0.53;
    %line([Vfinal Vfinal],[0 1]);
    %line([Vini Vini],[0 1]);
   
end    
%T = 1:1:length(M(:,1));
ax = gca;
xlim(ax,[0 1.2]);
xlabel('Binding Avidity (V)');
ylabel('Relative fitness at a given day');
%zlabel('Reproductive Number');

end