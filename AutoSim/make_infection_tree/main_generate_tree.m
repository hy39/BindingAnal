function [treeFile] = main_generate_tree(infile, smpno, starttime, endtime, display, lang)
%% Reconstruct phylogenies
% Without saving the files
% example1: main_generate_tree('voutput_small', '30', '10', '200', '0')
%
%% with antigenicTot
% example2: main_generate_tree(infile, '250', '30', '2000', '0', 'c')
% use example2 on PC based matlab
% produce phylogeny in nexus format and also a matlab file
%% 
% example4: lang = 's'
% quick and simple way to produce the tree
% won't save the matlab file
% this needs to be fixed

% Hsiang-Yu Yuan
% 12/01/2016

%%%p = path;
%%%p = path(p,'lib/');
%clear all; close all;
global epi_params;

if exist('infile','var') 
  if strcmp(infile,'')
      infile = ['dat/mod01/voutput_small'];
  else
  end
end

if isnumeric(smpno)
    n_seqs = smpno;
else
    n_seqs = str2num(smpno);
end

lg = 'c'; %c: c code; m:matlab; s:simple version in matlab 
epi_params.annotation='annotation'; % the 1st displayed annotation text
epi_params.display=0; % display the matlab figure?
epi_params.savefigure=1; % save the figure?
epi_params.savetree=2; % files to be saved. 0:no files, 1:only nexus tree file, 2:all the tree files
epi_params.lg=lg; % language platform

if exist('lang','var') 
  lg = lang;
end

if exist('display','var') 
  if isnumeric(display)
    epi_params.display=display;
  else
  if str2num(display) == 1
    epi_params.display=1;
  end
  end
end

if strcmp(lg,'c')
  epi_params.annotation='AntigenicDrift';
  %infile example: infile = ['dat/mod01/voutput_small']; 
  M = csvread([infile '.csv'],1);
  filename_infectionTreeData = strcat(infile, '_tree_', num2str(n_seqs), '_', num2str(starttime), '_', num2str(endtime));
  filename_infectionTreeData = strrep(filename_infectionTreeData,'input','output');
  dat_VirusesArray = M;
  count = length(dat_VirusesArray(:,1));%vid
  births = dat_VirusesArray(:,2);       %birth
  deaths = dat_VirusesArray(:,3);       %death
  parent = dat_VirusesArray(:,4);       %parentid
  infectionK = dat_VirusesArray(:,11);  %immnuity K
  binding = dat_VirusesArray(:,5);      %binding ini
  bindingFinal = dat_VirusesArray(:,6); %binding final
  antigenicTot = dat_VirusesArray(:,13);%total antigenic change
  %other improtant viral traits
  immJ = dat_VirusesArray(:,12);        %infected host immunity J
  parent(1:100) = 0;
end

if strcmp(lg,'m')
  %infile example: traitfile = ['dat/' proj '/virus_traits'];
  load(infile);
  filename_infectionTreeData = strcat(infile, '_tree_', num2str(n_seqs));
  filename_infectionTreeData = strrep(filename_infectionTreeData,'input','output');
  count = length(dat_VirusesArray(:,1));
  births = dat_VirusesArray(:,2);
  deaths = dat_VirusesArray(:,3);
  parent = dat_VirusesArray(:,4);
  infectionK = dat_VirusesArray(:,5);
  binding = dat_VirusesArray(:,7);
  bindingFinal = dat_VirusesArray(:,8); 
  antigenicTot = dat_VirusesArray(:,11); %%column numbers need to be changed
end

if strcmp(lg,'s')
  %[births, deaths, parent] = GetInfectionTree(infile);
  count = length(dat_VirusesArray(:,1));
  births = dat_VirusesArray(:,2);
  deaths = dat_VirusesArray(:,3);
  parent = dat_VirusesArray(:,4);  
end


if exist('starttime','var') 
    if isnumeric(starttime)
    % do nothing
    else
        starttime = str2num(starttime);
    end
else
    starttime = 30;
end
if exist('endtime','var')
    if isnumeric(endtime)
    % do nothing
    else
        endtime = str2num(endtime);
        if endtime < 1
            endtime = max(deaths)-1; 
        end
    end
else
    endtime = max(deaths)-1;
end
epi_params.tRange_stoch(1,1)=starttime;
epi_params.tRange_stoch(1,2)=endtime; %final version 365*45 (1968-2013)


locs = find(parent == 0); 
parent(locs) = NaN;

if strcmp(lg,'c')
    save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'infectionK', 'binding', 'bindingFinal','antigenicTot');
elseif strcmp(lg,'m')
    save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'infectionK', 'binding', 'bindingFinal');
end

treeFile = BuildTree_indiv_nexus(filename_infectionTreeData, n_seqs); % build nexus tree
% delete the temporary files
if exist([filename_infectionTreeData '.mat'], 'file') == 2
    %delete([filename_infectionTreeData '.mat']);
    %disp(['delete the temporary file: ' filename_infectionTreeData '.mat']);
end
end


