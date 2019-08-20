# SchoppikNagelLisberger2008
Code from Schoppik, Nagel &amp; Lisberger published in Neuron 2008

All relevant MATLAB code (92 KB) summarized briefly below. Functions not described are support functions

prepare.m --
Will re-create the unit???.mat files that comprise the analyzed data for each cell. 
Calls over 20GB of original files; if you really want to run this, please let me know 
and I'll make the original trial files available. I'd have to go dig through old hard drives, so <ugh>

aggregate.m -- extracts the relevant data from multiple single units, requires the unit???.mat files

aggregate_pairs_tbtv.m -- extracts the relevant multiple pairs, requires the pair???.mat files

smartCov.m -- generates correlation or covariation matrices with a known mean and variance

mp.m,mp2.m,sm.m -- will generate the figures from the output of aggregate.m

simulateResidualModel.m -- runs the model that comprises the supplementary figure
