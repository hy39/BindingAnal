% plot population fitness and population immunity
% Binding Scenario
% Fixed binding
% infile:
% infileV
% infile = 'dat/single_adaptive/.009/hostKs_1.csv'
% infileV = 'dat/single_adaptive/.009/voutput1_1.csv'

% infile = 'dat/single_adaptive_high/.004/hostKs_1.csv'
% infileV = 'dat/single_adaptive_high/.004/voutput1_1.csv'
function [filename_out] = plotPopRfromSkwithBinding(infile, infileV)

M = csvread( infile ); 

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
tii = 0;
totalN = sum(M(1,:));
T = [];
for t=1:1:length(M(:,1))
    %t = 81; %test
 

    tii = tii+1;
    T(tii) = t;
    %-- draw R curve only every 5 time units
    Sk = M(t,1:50);
    
    X = Sk;
    X = X./sum(X);
    %X = X./totalN;
    sk = X;
    f = @(v)-getOptV( v, del, sk );
    %Vopt(tii) = fminunc(f,v1,v2);
    Vopt(tii) = fminunc(f,[v2]);
    Vnaive(tii) = 0;
    Ropt(tii) = getPopR( del, sk, Vopt(tii));
    Rnaive(tii) = getPopR( del, sk, Vnaive(tii));

    
    if rem(t,10) == 1
    ti = ti + 1;
    
    

    
    %%----------------------------
    
    
    V = 0:0.02:1.5;
    for i=1:length(V)
        v = V(i);
        %R = getPopR( del, sk, vopt);
        R(i) = getPopR( del, sk, v);
    end
    
    plot3(V, -t*5*(repmat(1,length(V))), R, 'Color',[0.4,0.4,0.4]);
    end
end
%T = 1:1:length(M(:,1));
plot3(Vopt, -T*5, zeros(1,length(Vopt)), 'Color',[0.4,0.4,0.4]); % plot optimum V
plot3(1.5*(repmat(1,length(Vopt))), -T*5, Ropt, 'Color',[0.4,0.4,0.4]); % plot optimum R
plot3(1.5*(repmat(1,length(Vopt))), -T*5, Rnaive, 'Color',[0.4,0.4,0.4]); % plot naive R

%infileV = 'dat/voutput1_2';
%infileV = 'dat/slow/voutput1_2';
[Vini Vfinal Vmed] = calculateBinding( infileV, 5);
Vini = [Vini Vini(end)];
Vfinal = [Vfinal Vfinal(end)]; 
filename_out = infile(1:regexp(infile,'\.')-1);
filename_out = strrep(filename_out,'input','output');
filename_out = [filename_out '_binding.mat'];
save(filename_out,'Vfinal','Vini','Vopt');

Vmed = [Vmed Vmed(end)];
plot3(Vini, -T(1:length(Vini))*5, zeros(1,length(Vini)), 'R-'); % plot V inital
plot3(Vfinal, -T(1:length(Vini))*5, zeros(1,length(Vini)), 'B-'); % plot V final

ax = gca;
xlim
ax.YTick = [-5*length(M(:,1)):50:0];
ax.YTickLabel = [T(end)*5:-50:0];
ax.XMinorGrid = 'on'
ax.ZMinorGrid = 'on'
xlabel('Binding Avidity');
ylabel('Time (days)');
zlabel('Reproductive Number');

close(gcf);

end