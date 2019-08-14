
% Calculate p-values of Fisher's test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir = 'data\analysis\single_adaptive_high\005\';

load( strcat(dir,'NetChargeDiff.mat') );

threshold = 0.3; numCases = 200;

hasChange = sum((intMat<=threshold*-1),2) + sum((intMat>=threshold),2);
intTab = [sum(~isnan(intMat),2)-hasChange hasChange];
hasChange = sum((extMat<=threshold*-1),2) + sum((extMat>=threshold),2);
extTab = [sum(~isnan(extMat),2)-hasChange hasChange];
csvwrite('intTab.csv',intTab);
csvwrite('extTab.csv',extTab);

extTab = csvread('extTab.csv');
intTab = csvread('intTab.csv');

p = zeros(numCases,1);
for i = 1:numCases
    tab = [intTab(i,:); extTab(i,:)];
    [h,p(i),stats] = fishertest(tab);
end

csvwrite(strcat(dir,'p-values.csv'),p);