% plot the selection coefficient with mutational effects under different
% Binding Scenario
% Fixed binding
% example: infile = 'dat/lowV0/hostKs_low_adaptive.csv';

function [  ] = plotPopSCwithBinding(infile)

fixedV = 0.8;
MaxK = 11;
%vk = getVp(MaxK)

if exist('infile')
  M = csvread( infile );
  X1 = M(1,1:58)/sum(M(1,1:58));
  X1(1) = sum(X1(1:9));
  X1(2:9) = [];
  X2 = M(30,1:50)/sum(M(30,1:50)); %200 days
  X3 = M(200,1:50)/sum(M(200,1:50)); %300 days
else
X1 = normpdf(1:MaxK+1,1); % Weak Immunity       Sk(1) 
X2 = normpdf(1:MaxK+1,5); % Medium Immunity     sk(200) 
X3 = normpdf(1:MaxK+1,8); % Strong Immunity     Sk(300) 
end
delta = 1:1:15;
X = X1;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;

v1 = 0;
v2 = 1;
del = 0;
f = @(v)-getOptV( v, del, sk );
vopt = fminunc(f,v1,v2);



s = getPopSC( delta, sk, vopt)
%getOptV(delta,sk,v)
figure;
subplot(1,2,1);
plot(s);
hold on;

X = X2;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
v1 = 0; %calculate initial V
v2 = 1;
del = 0;
f = @(v)-getOptV( v, del, sk );
vopt = fminunc(f,v1,v2);
s = getPopSC( delta, sk, vopt)
plot(s);

X = X3;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
v1 = 0; %calculate initial V
v2 = 1;
del = 0;
f = @(v)-getOptV( v, del, sk );
vopt = fminunc(f,v1,v2);
s = getPopSC( delta, sk, vopt)
plot(s);

%-------------------------------------------------------------------------
subplot(1,2,2);
X = X1;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
new_delta = [0 delta];
for d=1:length(new_delta)
    v1 = 0;
    v2 = 1;
    f = @(v)-getOptV( v, new_delta(d), sk );
    del = [0 new_delta(d)];
    vopt = fminunc(f,v1,v2);
    s_tmp = getPopSC( del, sk, vopt-0.2);
    s(d) = s_tmp(2);
end
plot(s);
hold on;

X = X2;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
for d=1:length(new_delta)
    v1 = 0;
    v2 = 1;
    f = @(v)-getOptV( v, new_delta(d), sk );
    del = [0 new_delta(d)];
    vopt = fminunc(f,v1,v2);
    s_tmp = getPopSC( del, sk, vopt-0.2);
    s(d) = s_tmp(2);
end
plot(s);

X = X3;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
for d=1:length(new_delta)
    v1 = 0;
    v2 = 1;
    f = @(v)-getOptV( v, new_delta(d), sk );
    del = [0 new_delta(d)];
    vopt = fminunc(f,v1,v2);
    s_tmp = getPopSC( del, sk, vopt-0.2);
    s(d) = s_tmp(2);
end
plot(s);
xlim = [0 15]





end