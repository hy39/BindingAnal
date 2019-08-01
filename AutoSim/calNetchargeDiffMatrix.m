


%% II.

% infileV1 - e.g. data\output\single_adaptive_high\004\voutput1_1_tree_300_50_415_simple.mat
% infileH1 - e.g. data\output\single_adaptive_high\004\hostKs_1_binding.mat
function[] = calNetchargeDiffMatrix(numCases, smpno, infileV1, infileH1)
intMat = zeros(1,smpno); extMat = zeros(1,smpno);
startIndex = 1;
marker = '_1_';
for i = startIndex:numCases
    replace = strcat('_',num2str(i),'_');
    infileV = strrep(infileV1, marker, replace);
    infileH = strrep(infileH1, marker, replace); 
    [intMat, extMat] = netChargeByTimeSim_2out(infileV, infileH, intMat, extMat);
end
intMat(1,:) = []; extMat(1,:) = [];

% configure savefile path and then save
saveDir = extractBefore(infileH1, max( strfind(infileH1,'\') ) ); % returns the directory of infileH1
saveDir = strrep(saveDir,'output','analysis'); saveDir = strrep(saveDir,'input','analysis');
saveFile = strcat(saveDir,'\NetChargeDiff.mat');
save(saveFile, 'intMat', 'extMat');
end