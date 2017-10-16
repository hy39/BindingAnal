20150526
The tradeoff package demonstrate the withinhost R0 and transmission rate given different immune status, virus binding avidity
and other parameters values. 

Check the immunehistory_model.pdf document in antigenic_drift/doc/tradeoff for more details.

To plot reproductive number and the optimum binding avidity:
workdir E:\Documents\Github\BindingAnal\fitness\tradeoff\Mtradeoff
function: plotPopRfromSkwithBinding(infile, infileV)
example:
infileV = 'dat/lowV0/voutput1_low_adaptive.csv';
infile = 'dat/lowV0/hostKs_low_adaptive.csv'
plotPopRfromSkwithBinding(infile, infileV), this file plot the population level reproductive number R
