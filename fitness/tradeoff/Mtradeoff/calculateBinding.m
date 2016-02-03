function [ averageVini averageVfinal] = calculateBinding( infileV, interval )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  M = csvread([infileV '.csv'],1);
  dat_VirusesArray = M;
  count = length(dat_VirusesArray(:,1));%vid
  births = dat_VirusesArray(:,2);       %birth
  deaths = dat_VirusesArray(:,3);       %death
  parent = dat_VirusesArray(:,4);       %parentid
  infectionK = dat_VirusesArray(:,11);  %immJ
  binding = dat_VirusesArray(:,5);      %binding ini
  bindingFinal = dat_VirusesArray(:,6); %binding final
  antigenicTot = dat_VirusesArray(:,12);%total antigenic change
  averageVini = 0;
  averageVfinal = 0;
  i=0;
  for t=1:interval:300
    i = i + 1
    averageVini(i) = mean(binding(find(births<t & deaths>t)));
    averageVfinal(i) = mean(bindingFinal(find(births<t & deaths>t)));
  end
end
