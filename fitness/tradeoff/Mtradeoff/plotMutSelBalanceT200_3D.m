% plot the selection coefficient with mutational effects under different
% Binding Scenario
% Fixed binding
% Ploat the fitness surface with different antigenic changes
% The function was replaced by plotMutSelBalance_3D()

function [  ] = plotMutSelBalanceT200_3D()
%infile = 'dat/lowV0/hostKs_low_adaptive.csv';
infile = 'dat/single_adaptive_high/.004/hostKs_1.csv';
M = csvread( infile ); 
%V = csvread('dat/lowV0/voutput1_low_adaptive.csv',2);

delta = 1:15;


del = 10; %antigenic change 2 as our simulation


%figure;
hold on;
ti=0;
tii = 0;
totalN = sum(M(1,:));
T = [];
j = 1;
for t=1:10:length(M(:,1))

    %-- draw R curve only every 5 time units
    %test 1 200 
    
    
    t = 200;
    Sk = M(t,1:50);
    
    X = Sk;
    %X = X./sum(X);
    X = X./totalN;
    sk = X;
    f = @(v)-getOptV( v, del, sk );
    v1 = 0.2;
    v2 = 1;
    %Vopt = fminunc(f,[v1,v2]);
    Vopt = fminunc(f,[v2]);
    Vnaive = 0;
    Ropt = getPopR( del, sk, Vopt);
    Rnaive = getPopR( del, sk, Vnaive);

    %Calculate current Virus
    
    %%----------------------------

    V = 0.1:0.05:1.5;
    
    
    for i=1:length(V)
       %for j=1:1:10
        %del = j-1;
        v = V(i);
        %R = getPopR( del, sk, vopt);
        Rpop(j,i) = getPopR( del, sk, v);
        %Rin(i) = getPopRin( del, sk, v);
        Inewt(i,:) = getBeta( sk, del, v ).*sk; % estimate the newly infected individuals
        Rin(j,i) = getPopRin( del, Inewt(i,:), v);
       %end
    end
    j = j + 1;
end

    x = 0.1:0.05:1.5;
    y = 1:10:length(M(:,1));
    [X,Y] = meshgrid(x,y);
    surf(X,Y,Rpop);
    hold on;
    surf(X,Y,Rin);
    %Beta = R*(1/3.3);
    
    
    %plot3(V, -t*5*(repmat(1,length(V))), R, 'Color',[0.4,0.4,0.4]);
    %plot3(2*ones(1,length(V)),V, Rin/max(Rin));
    %hold on;
    %plot3(2*ones(1,length(V)),V, R/max(R));
    %plot(V, Rin/max(Rin));
    %plot(V, R/max(Rpop));

  
ax = gca;
xlim(ax,[0 1.2]);
xlabel('Binding Avidity (V)');
ylabel('Relative fitness at a given day');
%zlabel('Reproductive Number');

end