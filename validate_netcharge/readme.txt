20160127
Plot the netcharge vs age
1. go to the folder BindingAnal/NetCharge/mat/Mnetcharge/dat/h3n2_ny run genTable and produce Chlist.dat
  - Binomial confident interval is calculated
2. Back to the currect folder and run netcharge_age.R


20160122
Plot the correlation of binding and netcharge
1. save validate_netcharge.xlsl to
	a. RDE/dat_netch_binding_hensley.csv
	b. RDE/dat_netch_binding_others.csv
	c. RDE/dat_netch_binding_total.csv
2. plotnetbinding.m transform dat_netch_binding_hensley.csv to MHensley.dat
3. plot_bindingchanges.R read MHensley.dat and plot the ressult


20160119
plotbetbinding.m would save the netcharge and binding data as csv file
plot_bindingchanges.R would plot the distribution of binding versus netcharge. 


what about h1n1?