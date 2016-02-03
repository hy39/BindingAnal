% plot the selection coefficient with mutational effects under different
% Binding Scenario
% Fixed binding
% 

function [  ] = plotPopSCwithBinding( )

fixedV = 0.8;
MaxK = 11;
%vk = getVp(MaxK)

X1 = normpdf(1:MaxK+1,1);
X2 = normpdf(1:MaxK+1,5);
X3 = normpdf(1:MaxK+1,8);
delta = 1:15;
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
    s_tmp = getPopSC( del, sk, vopt);
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
    s_tmp = getPopSC( del, sk, vopt);
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
    s_tmp = getPopSC( del, sk, vopt);
    s(d) = s_tmp(2);
end
plot(s);
xlim = [0 15]





end