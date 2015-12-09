function [ output_args ] = main_plot_figure1( age_min, age_max, saveoutput )
%PLOT_BINDING_REGRESSION Summary of this function goes here
%   Detailed explanation goes here

  if ~exist('saveoutput')
      saveoutput = '';
  end
  hFig = figure;
  Figw = 380*1.3;
  Figh = 630*1.9;
  set(hFig, 'Position', [100 100 Figw Figh]);
  subplot(3,1,1);
  plot_net_byage_reg(age_min, age_max, saveoutput);
  
  subplot(3,1,2);
  plot_bding_byage_reg(age_min, age_max, saveoutput)
 
  subplot(3,1,3);
  plot_h1n1_net_byage_reg(age_min, age_max, saveoutput)
end

