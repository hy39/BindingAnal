% plot the selection coefficient with mutational effects under different
% Binding Scenario
% Fixed binding
% 

function [  ] = plotPopRfromSkwithBinding(infile)

M = csvread( infile ) 

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
del = 0;
%f = @(v)-getOptV( v, del, sk );
%vopt = fminunc(f,v1,v2);

figure;
hold on;
ti=0;

for t=1:1:length(M(:,1))
    Sk = M(t,1:50);
    X = Sk;
    X = X./sum(X);
    sk = X;
    f = @(v)-getOptV( v, del, sk );
    Vopt(t) = fminunc(f,v1,v2);
    Vnaive(t) = 0;
    Ropt(t) = getPopR( del, sk, Vopt(t));
    Rnaive(t) = getPopR( del, sk, Vnaive(t));
 
    %-- draw R curve only every 5 time units
    if rem(t,5) == 1
    ti = ti + 1;
    V = 0:0.02:1.5;
    for i=1:length(V)
        v = V(i);
        %R = getPopR( del, sk, vopt);
        R(i) = getPopR( del, sk, v);
    end
    
    plot3(V, -t*5*(repmat(1,length(V))), R, 'Color',[0.4,0.4,0.4]);
    end
end
T = 1:1:length(M(:,1));
plot3(Vopt, -T*5, zeros(1,length(Vopt)), 'Color',[0.4,0.4,0.4]); % plot optimum V
plot3(1.5*(repmat(1,length(Vopt))), -T*5, Ropt, 'Color',[0.4,0.4,0.4]); % plot optimum R
plot3(1.5*(repmat(1,length(Vopt))), -T*5, Rnaive, 'Color',[0.4,0.4,0.4]); % plot naive R

infileV = 'dat/voutput1_2';
[Vini Vfinal] = calculateBinding( infileV, 5);
 
plot3(Vini, -T*5, zeros(1,length(Vopt)), 'R-'); % plot V inital
plot3(Vfinal, -T*5, zeros(1,length(Vopt)), 'B-'); % plot V final

ax = gca;
xlim
ax.YTick = [-350:50:0];
ax.YTickLabel = [350:-50:0];
ax.XMinorGrid = 'on'
ax.ZMinorGrid = 'on'
xlabel('Binding Avidity');
ylabel('Time (days)');
zlabel('Reproductive Number');

end