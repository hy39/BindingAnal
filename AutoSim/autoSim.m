


% - for now, input filenames, e.g. voutput1, hostKs, etc. will be hard-coded
% - for now, filenames' index x will also be hard-coded, e.g.
% voutput1_x.csv, hostKs_x.csv

% I. generate all outputs - phylogeny tree, host binding, binding avidity
% curves
% - input parameters: # of cases, # of samples for each case, starttime for
% each case, endtime for each case

% II. further analysis - avgBinding, calNetchargeMatrix

% III. plots

% XX. further analysis - Fisher's tests' p-values, selection, cumulative
% incidence



%% I. 

% infileVHeader - e.g. data\input\single_adaptive_high\004\voutput1 
% infileHHeader - e.g. data\input\single_adaptive_high\004\hostKs
% (NOTE: use '\' instead of '/'!)
% (NOTE: if the input files are copied directly from the output of the 'bindingavid'
% project at R, then do not modify their names, i.e. voutput1 and hostKs.
function[] = autoSim(numCases, infileVHeader, smpno, starttime, endtime, infileHHeader)
startIndex = 1;
for i = startIndex:numCases
    infileV = strcat(infileVHeader, '_', num2str(i));
    treeFile = main_generate_tree(infileV, smpno, starttime, endtime, '0', 'c');
    infileH = strcat(infileHHeader, '_', num2str(i));
    bindingFile = plotPopRfromSkwithBinding( strcat(infileH,'.csv') , strcat(infileV,'.csv') );
    plotNetchargeByTimeSim(treeFile, bindingFile);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


