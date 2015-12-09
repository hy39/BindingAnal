function [] = datafitting( )

load('datafitting_clean.mat');
virus = virus8_11;
%virus = virus(find(virus(:,1)>0.5),:)
figure;
X = virus(:,1)
Y = virus(:,2)
x = linspace(min(X),max(X));
scatter(X,Y,'k')
hold on;
span = 0.1;
%line(x,mylowess([X,Y],x,span))


scatter(X, Y)
f = @(xy) mylowess(xy,x,span);
yboot2 = bootstrp(500,f,[X,Y])';
meanloess = mean(yboot2,2);
h1 = line(x, meanloess,'color','k','linestyle','-','linewidth',2);

stdloess = std(yboot2,0,2);
h2 = line(x, meanloess+2*stdloess,'color','r','linestyle','--','linewidth',2);
h3 = line(x, meanloess-2*stdloess,'color','r','linestyle','--','linewidth',2);

L5 = legend('Localized Regression','Confidence Intervals',2);
snapnow
end

