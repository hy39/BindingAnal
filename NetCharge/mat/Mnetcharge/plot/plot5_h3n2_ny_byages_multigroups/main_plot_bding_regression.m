function [ output_args ] = main_plot_bding_regression( age_min, age_max, saveoutput )
%PLOT_BINDING_REGRESSION Summary of this function goes here
%   Detailed explanation goes here

  if ~exist('saveoutput')
      saveoutput = '';
  end
  hFig = figure;
  Figw = 380;
  Figh = 630;
  set(hFig, 'Position', [100 100 Figw Figh]);
  subplot(2,1,1);
  plot_net_byage_reg(age_min, age_max, saveoutput);
  
  subplot(2,1,2);
  plot_bding_byage_reg(age_min, age_max, saveoutput)
  
  %subplot(3,1,3);
  %plot_h1n1_net_byage_reg(age_min, age_max, saveoutput)
  
  %subplot(2,2,4);
  %plot_h1n1_bding_byage_reg (age_min, age_max, saveoutput)
end

