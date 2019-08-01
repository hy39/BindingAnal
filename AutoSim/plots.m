
dir = 'data\analysis\single_adaptive_high_custom2\012\';



%% III.

load(strcat(dir,'VMat.mat'));
VintMat(VintMat==0) = NaN;
VextMat(VextMat==0) = NaN;
avgVint = nanmean(VintMat(:,2:end),2);
avgVext = nanmean(VextMat(:,2:end),2);
avgVopt = nanmean(VoptMat(:,2:end),2);
stdVint = std(VintMat(:,2:end)','omitnan')';
stdVext = std(VextMat(:,2:end)','omitnan')';
stdVopt = std(VoptMat(:,2:end)','omitnan')';

% mean 
figure
hold on;
plot(VintMat(:,1),avgVint,'b.' , VextMat(:,1),avgVext,'gx');
plot(VoptMat(:,1),avgVopt,'-','Color',[0.4 0.4 0.4]);
ylabel('Binding Avidity (V)');
xlabel('Time (days)');
legend({'internal','external'}, 'location','southeast');
hold off;

% mean & 95% CI
figure
hold on;
plot(VintMat(:,1),avgVint,'b' , VextMat(:,1),avgVext,'g');
plot(VoptMat(:,1),avgVopt,'k-');
fill([VintMat(2:end,1); flipud(VintMat(2:end,1))], [avgVint(2:end)-1.96.*stdVint(2:end); flipud(avgVint(2:end)+1.96.*stdVint(2:end))], 'b', 'facecolor','blue', 'edgecolor','none', 'facealpha','0.3'); 
fill([VextMat(4:end,1); flipud(VextMat(4:end,1))], [avgVext(4:end)-1.96.*stdVext(4:end); flipud(avgVext(4:end)+1.96.*stdVext(4:end))], 'g', 'facecolor','green', 'edgecolor','none', 'facealpha','0.3'); 
fill([VoptMat(:,1); flipud(VoptMat(:,1))], [avgVopt-1.96.*stdVopt; flipud(avgVopt+1.96.*stdVopt)], 'k', 'facecolor',[0.4 0.4 0.4], 'edgecolor','none', 'facealpha','0.3'); 
ylabel('Binding Avidity (V)');
xlabel('Time (days)');
axis([0 415 0 1]);
legend({'internal','external','optimum'}, 'location','southeast');
hold off;

% individual points
figure
hold on;
for i=2:size(VintMat,1)
   plot(VintMat(:,1),VintMat(:,i),'b.');
   plot(VextMat(:,1),VextMat(:,i),'gx');
end
for i=2:size(VoptMat,1)
    plot(VoptMat(:,1),VoptMat(:,i),'-','Color',[0.4 0.4 0.4]);
end
ylabel('Binding Avidity (V)');
xlabel('Time (days)');
legend({'internal','external'}, 'location','southeast');
hold off;



%% IV.

load(strcat(dir,'NetChargeDiff.mat'));
edges = -1:0.05:1;
centres = edges(1:end-1)+ diff(edges)/2;
intBin = zeros(1,length(edges)-1); extBin = zeros(1,length(edges)-1);
for i=1:size(intMat,1)
   intBin = intBin + histcounts(intMat(i,:), edges);
   extBin = extBin + histcounts(extMat(i,:), edges);
end
figure
hold on;
plot(centres, intBin, 'b');
plot(centres, extBin, 'g');
xlabel('Charge difference from ancestor');
ylabel('No. of nodes ');
legend({'internal','external'}, 'location','northeast');
hold off;
