


%% II.

% infileV1 - e.g. data\output\single_adaptive_high\004\voutput1_1_tree_300_50_415_simple.mat
% infileH1 - e.g. data\output\single_adaptive_high\004\hostKs_1_binding.mat
function[] = avgBinding(numCases, infileV1, infileH1)
startIndex = 1;
VintMat = [10:15:545]'; % NOTE: the ending time corresponds to the argument endtime for autoSim
VextMat = VintMat;
VoptMat = [10:5:545]';
marker = '_1_';
for i = startIndex:numCases
    replace = strcat('_',num2str(i),'_');
    infileV = strrep(infileV1, marker, replace);
    infileH = strrep(infileH1, marker, replace);
    [VintMat, VextMat, VoptMat] = netChargeByTimeSim_3out(infileV, infileH, VintMat, VextMat, VoptMat);   
end

% configure savefile path and then save
saveDir = extractBefore(infileH1, max( strfind(infileH1,'\') ) ); % returns the directory of infileH1
saveDir = strrep(saveDir,'output','analysis'); saveDir = strrep(saveDir,'input','analysis');
saveFile = strcat(saveDir,'\VMat.mat');
save(saveFile, 'VintMat', 'VextMat', 'VoptMat');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

