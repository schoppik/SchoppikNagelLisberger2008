%
%
%
%
%          by: david schoppik
%        date: 08-15-2006
%     purpose: To create a population data structure ("data") for
%              statistical analysis from the "unit???.mat" files
%

tic

% variables
data_dir = '/home/schoppik/exp/analyzeddata/';
data = [];
smfun = normalizedgaussian([0 20],-100:100);
prefdirvec = 45:45:360;
aboveflankvec = circshift(45:45:360,[0 -1]);
belowflankvec = circshift(45:45:360,[0 1]);

% pulls the ??? from the unit???.mat files; assumes a certain number of
% non-unit???.mat files; 

cd(data_dir)
d = dir;
d = vertcat(d(5:end).name);
unit_list = str2num(d(:,5:7));

for unit_ind = 1:length(unit_list);
  disp(unit_ind)
  
  % load the unit and place it into a placeholder
  load(sprintf('unit%.3d',unit_list(unit_ind)));
  eval(sprintf('u = unit%.3d;',unit_list(unit_ind)));
  
  %%%%%%%%%%%%
  % Raw data %
  %%%%%%%%%%%%
  
  data.binaryspikes{unit_ind} = u.binaryspikes;
  data.evel{unit_ind} = u.evel;
  
  %%%%%%%%%%%
  % Filters %
  %%%%%%%%%%%
  
  % get the filter
  data.filters.rev(unit_ind,:) = mean(u.filterbank.var,2);
  data.filters.fwd(unit_ind,:) = mean(u.forward.filterbank.var,2);
  data.filters.abv(unit_ind,:) = mean(u.above.filterbank.var,2);
  data.filters.blw(unit_ind,:) = mean(u.below.filterbank.var,2);
    
  % evaluate the filter
  nfft = 512;
  Frand = abs(fft(u.filterbank.randvar.^2,nfft));
  Ffilt = abs(fft(u.filterbank.var.^2,nfft));
  
  data.filters.randpw(unit_ind,:) = sum(Frand);
  data.filters.filtpw(unit_ind,:) = sum(Ffilt);

  % how do the forward and reverse filters compare?
  % get the average "flippability" of each filter
  Fr = data.filters.rev(unit_ind,:)/norm(data.filters.rev(unit_ind,:));
  Ff = data.filters.fwd(unit_ind,:)/norm(data.filters.fwd(unit_ind,:));
  
  data.filters.similar(unit_ind) = Fr*fliplr(Ff)';

  % compute the coherence
  nfft = 2048;
  A = u.evel - repmat(mean(u.evel),[size(u.evel,1) 1]);
  B = u.binaryspikes - repmat(mean(u.binaryspikes),[size(u.binaryspikes,1) 1]);
  fAB = mean(fft(A',nfft).*conj(fft(B',nfft)),2);
  fAA = mean(fft(A',nfft).*conj(fft(A',nfft)),2);
  fBB = mean(fft(B',nfft).*conj(fft(B',nfft)),2);
  data.coherence(unit_ind,:) = (abs(fAB).^2)./(fAA.*fBB);
  
  % how well do the filters generalize to data from another direction?
  % we'll express it as a ratio, so let's do it in the preferred direction
  % first
  
  filt = mean(u.filterbank.var,2);
  dex = find(u.alltype == prefdirvec(find(prefdirvec == u.prefdir)));
  spvar = u.allspikes(dex,:) - repmat(mean(u.allspikes(dex,:)),[length(dex) 1]);
  evar = u.allevel(dex,:) - repmat(mean(u.allevel(dex,:)),[length(dex) 1]);
  ehat = myconv(spvar',repmat(filt,[1 length(dex)]));
  
  % first, using the complete way of assesing the data
  data.crosspred.pd(unit_ind) = mean(mycorr(evar',ehat));
  
  % now, using Steve's maximum
  data.crosspred.maxpd(unit_ind) = max(mycorr(evar,ehat'));
  
  % now, we've got to be careful to flip the filter appropriately.
  % if it was originally generated on 180 or 270, invert it with "-"
  % because the next two tests will be on flanks that are positive.
  
  if u.prefdir == 180 | u.prefdir == 270
    filt = -filt;
  end
  
  dex = find(u.alltype == aboveflankvec(find(prefdirvec == u.prefdir)));
  spvar = u.allspikes(dex,:) - repmat(mean(u.allspikes(dex,:)),[length(dex) 1]);
  evar = u.allevel(dex,:) - repmat(mean(u.allevel(dex,:)),[length(dex) 1]);
  
  % if we're testing on 180 or 270, we've got to invert the filter as well
  if aboveflankvec(find(prefdirvec == u.prefdir)) == 180 |...
      aboveflankvec(find(prefdirvec == u.prefdir)) == 270
    filt = -filt;
  end
  
  ehat = myconv(spvar',repmat(filt,[1 length(dex)]));
  
  % first, using the complete metric
  data.crosspred.above(unit_ind) = mean(mycorr(evar',ehat));
  
  % now Steve's
  data.crosspred.maxabove(unit_ind) = max(mycorr(evar,ehat'));

  % now check how the filter that was generated on that data would do on
  % the same metric
  ehat = myconv(spvar',repmat(mean(u.above.filterbank.var,2),[1 length(dex)]));
  data.crosspred.actabove(unit_ind) = max(mycorr(evar,ehat'));

  
  % if we're testing on 180 or 270, we've got to flip it back, otherwise
  % the next bit of testing will be messed up.  sloppy sloppy code.  sorry.
  if aboveflankvec(find(prefdirvec == u.prefdir)) == 180 |...
      aboveflankvec(find(prefdirvec == u.prefdir)) == 270
    filt = -filt;
  end
    
  
  % test on the below flank as well
  
  dex = find(u.alltype == belowflankvec(find(prefdirvec == u.prefdir)));
  spvar = u.allspikes(dex,:) - repmat(mean(u.allspikes(dex,:)),[length(dex) 1]);
  evar = u.allevel(dex,:) - repmat(mean(u.allevel(dex,:)),[length(dex) 1]);
   % if we're testing on 180 or 270, we've got to invert the filter as well
  if belowflankvec(find(prefdirvec == u.prefdir)) == 180 |...
      belowflankvec(find(prefdirvec == u.prefdir)) == 270
    filt = -filt;
  end
  
  ehat = myconv(spvar',repmat(filt,[1 length(dex)]));
 
  % first, the complete
  data.crosspred.below(unit_ind) = mean(mycorr(evar',ehat));
 
  % then, Steve's
  data.crosspred.maxbelow(unit_ind) = max(mycorr(evar,ehat'));
  
  % now check how the filter that was generated on that data would do on
  % the same metric
  ehat = myconv(spvar',repmat(mean(u.below.filterbank.var,2),[1 length(dex)]));
  data.crosspred.actbelow(unit_ind) = max(mycorr(evar,ehat'));
  
  %%%%%%%
  % REV %
  %%%%%%%
  
  % what was the corr?
  data.corr.rev(unit_ind) = mean(u.test.corr);
  data.shufcorr.rev(unit_ind) = mean(u.test.shufcorr);
  
  % was that significant?
  [data.corrsig.rev(unit_ind) data.corrp.rev(unit_ind)] = ttest(u.test.corr,u.test.shufcorr,0.5/(150*701));
  
  % how about as a function of time
  data.timecorr.rev(unit_ind,:) = mean(u.test.timecorr);
  data.shuftimecorr.rev(unit_ind,:) = mean(u.test.shuftimecorr);
  %data.timecorrsig.rev(unit_ind,:) = ttest(u.test.timecorr,u.test.shuftimecorr,0.05/(150*701),'right');
  %[data.timecorrsig.rev(unit_ind,:) data.timecorrp.rev(unit_ind,:)] = ttest(u.test.timecorr,u.test.shuftimecorr,0.5/(150*701),'right');
  data.timecorrsig.rev(unit_ind,:) = (mean(u.test.timecorr) > 1.5*std(u.test.shuftimecorr));
  
  
  % how did the filter do at predicting the mean?
  data.crosscorr.rev(unit_ind) = mean(u.cross.corr);
  data.bestcrosscorr.rev(unit_ind) = mean(u.cross.bestcorr);
  data.shufcrosscorr.rev(unit_ind) = mean(u.cross.shufcorr);
  
  % was the performance at predicting the mean significant?
  [data.crosssig.rev(unit_ind) data.crossp.rev(unit_ind)] = ttest(u.cross.corr,u.cross.shufcorr);
    
  %%%%%%%
  % FWD %
  %%%%%%%
  
  % now, the forward filters
  data.corr.fwd(unit_ind) = mean(u.forward.test.corr);
  data.shufcorr.fwd(unit_ind) = mean(u.forward.test.shufcorr);
  
  % was that significant?
  [data.corrsig.fwd(unit_ind) data.corrp.fwd(unit_ind)] = ttest(u.forward.test.corr,u.forward.test.shufcorr,0.5/(150*701));
  
  % how about as a function of time
  data.timecorr.fwd(unit_ind,:) = mean(u.forward.test.timecorr);
  data.shuftimecorr.fwd(unit_ind,:) = mean(u.forward.test.shuftimecorr);
  [data.timecorrsig.fwd(unit_ind,:) data.timecorrp.fwd(unit_ind,:)] = ttest(u.forward.test.timecorr,u.forward.test.shuftimecorr);
  
  % how did the filter do at predicting the mean?
  data.crosscorr.fwd(unit_ind) = mean(u.forward.cross.corr);
  data.bestcrosscorr.fwd(unit_ind) = mean(u.forward.cross.bestcorr);
  data.shufcrosscorr.fwd(unit_ind) = mean(u.forward.cross.shufcorr);
  
  % was that prediciton significant?
  [data.crosssig.fwd(unit_ind) data.crossp.fwd(unit_ind)] = ttest(u.forward.cross.corr,u.forward.cross.shufcorr);

  %%%%%%%%%
  % Other %
  %%%%%%%%%
  
  % aligned
  data.corr.aln(unit_ind) = mean(u.aligned.test.corr);
  data.timecorr.aln{unit_ind} = mean(u.aligned.test.timecorr);
  
  % short
  data.corr.srt(unit_ind) = mean(u.short.test.corr);
  data.timecorr.srt{unit_ind} = mean(u.short.test.timecorr);
  
  % above
  data.corr.abv(unit_ind) = mean(u.above.test.corr);
  data.timecorr.abv(unit_ind,:) = mean(u.above.test.timecorr);
  
  % below
  data.corr.blw(unit_ind) = mean(u.below.test.corr);
  data.timecorr.blw(unit_ind,:) = mean(u.below.test.timecorr);
  
  % position
  data.corr.pos(unit_ind) = mean(u.position.test.corr);
  data.timecorr.pos(unit_ind,:) = mean(u.position.test.timecorr);
  
  %%%%%%%%%%%%%%%%%%%
  % covariance data %
  %%%%%%%%%%%%%%%%%%%
  
  data.cov(unit_ind,:) = mean(u.cov.test.timecorr);
  
  %%%%%%%%%%%%%%%%
  % General Info %
  %%%%%%%%%%%%%%%%
   
  % numtrials
  data.numtrials(unit_ind) = size(u.binaryspikes,1);
  
  % preferred direction
  data.prefdir(unit_ind) = u.prefdir;
  
  % average firing rate
  data.fr(unit_ind) = 1000.*(sum(sum(u.binaryspikes(:,201:end))/prod(size(u.binaryspikes(:,201:end)))));

  % what was the unitid
  data.unitid(unit_ind) = unit_list(unit_ind);
  
  % what was the calculated preferred direction
  data.pd(unit_ind,:) = u.pd;
  
  % what was the PSTH
  data.psth(unit_ind,:) = sum(u.binaryspikes)./size(u.evel,1);
  
  frate = conv(data.psth(unit_ind,:),smfun);
  data.smoothpsth(unit_ind,:) = frate(101:end-100);
  
  % how much did the PSTH look like the performance of the model in time?
  data.timevpsth(unit_ind,:) = corr2(mean(u.test.timecorr),data.psth(unit_ind,:));
  
  % what was the total variance of the eye movement?  of the spike trains?
  %data.var.eye(unit_ind) = mean(var(u.evel'));
  %data.var.sp(unit_ind) = mean(var(u.binaryspikes'));
 
  data.var.sp(unit_ind) = var(sum(u.binaryspikes,2));
  
  % how well isolated was the unit?
  % data.rank{unit_ind} = u.rank;
  
  % how many spikes were there in each direction?
   data.spikesbydir(unit_ind,:) = u.sp;

  % clear out the unit's data (for memory reasons)
  clear(sprintf('unit%.3d',unit_list(unit_ind)));
end

toc


%%%%%%%%%%%%%%%%%%%%
% Summary Analyses %
%%%%%%%%%%%%%%%%%%%%

% test for significance
h = ttest(data.filters.randpw',data.filters.filtpw',.05/150);

data.sigind = find(h);
data.nsigind = find(~h);

% tweak the pd values so that they all fall to the right of 0, and are
% b/w 0 and 360 (sometimes the fitting algorithm converges b/w -360 and 0)
data.pd(find(data.pd(:,1) > 720),:) = [nan nan];
data.pd(find(data.pd(:,1) < -720),:) = [nan nan];
dex = find(data.pd(:,1) > 360);
data.pd(dex,1) = data.pd(dex,1) - 360;
dex = find(data.pd(:,1) < 0);
data.pd(dex,1) = data.pd(dex,1) + 360;

% the tuning bandwidth should be constrained to be 360 if it is larger
% i.e. 360 = untuned
data.pd(find(data.pd(:,2) > 360),2) = 360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter Characterization Analyses %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the next set of data was gotten by eye:
% primarily single-lobed
data.filters.type([3 19 20 24 26 37 40 51 56 60 62 77 90 96 97 102 109 115 124 139]) = 1;
% primarily double-lobed
data.filters.type([1 6 7 8 11 12 16 17 25 27 32 34 41 47 49 52 53 54 57:59 61 63 65 66 68 69,...
  71 72 82 89 100 104 110 111 118 128 137 138 140]) = 2;
% primarily triple-lobed
data.filters.type([2 4 5 18 38 87 97 107]) = 3;

% unclassifiable
data.filters.type([9 10 13 14 15 21:23 28:31 33 35 36 39 42:46 48 50 55 64 67 70 73:76 78:81 83:86,...
  88 91:95 98 99 103 105 106 109 101 108 112:114 116 117 119 120:123 125:127 129:136 141]) = nan;

% noteworthy (to me) = [63 87 112 118];

data.filters.lat = [-4 -13 -85 -57 -11,... % 5
  -56 79 -95 nan nan,... % 10
  11 2 nan nan nan,... % 15
  -10 -131 -110 -72 70,... % 20
  nan nan nan -17 -40,... % 25 
  -14 -28 nan nan nan,... % 30
  nan 16 nan -16 nan,... % 35
  nan -90 41 nan 25,... % 40
  -5 nan nan nan nan,... % 45
  nan 15 nan 150 nan,... % 50 
  5 -48 -111 -10 nan,... % 55
  9 -45 -30 8 55,... % 60
  8 -80 0 nan -33,... % 65
  -19 nan 10 17 nan,... % 70
  -28 -70 nan nan nan,... % 75
  nan -15 nan nan nan,... % 80
  nan -26 nan nan nan,... % 85
  nan 154 nan -30 34,... % 90
  nan nan nan nan nan,... % 95
  29 -140 nan nan -20,... % 100
  nan 30 nan -21 nan,... % 105
  nan -5 nan nan 3,... % 110
  -2 nan nan nan 10,... % 115
  nan nan 10 nan nan,... % 120
  nan nan nan 44 nan,... % 125
  nan nan 60 nan nan,... % 130 
  nan nan nan nan nan,... % 135
  nan 32 18 50 -4 nan]; % 141

% unclassifiable indices, for plotting
unk = find(isnan(data.filters.type));
data.unkind = unk(~ismember(unk,data.nsigind));
data.classsigind = data.sigind(~ismember(data.sigind,unk));

% Steve's preferred way of quantifying the neuron-behavior correlation
cutoff = 0.9;
f = [0:2047]*1000/2048;


for i = 1:141,
  data.summary.maxcorr.rev(i) = max(data.timecorr.rev(i,:));
  data.summary.maxcorr.blw(i) = max(data.timecorr.blw(i,:));
  data.summary.maxcorr.abv(i) = max(data.timecorr.abv(i,:));
  data.summary.maxcorr.pos(i) = max(data.timecorr.pos(i,:));
  data.summary.maxcorr.aln(i) = max(data.timecorr.aln{i});
  data.summary.maxcorr.srt(i) = max(data.timecorr.srt{i});
  
  pottime = diff(find(data.timecorr.rev(i,:) >= (data.summary.maxcorr.rev(i)*cutoff)));
  data.summary.maxcorrtime.rev(i) = length(find(pottime == 1));
  data.summary.maxcoh.rev(i) = max(data.coherence(i,1:100));
  data.summary.maxfreq.rev{i} = f(find(data.coherence(i,1:100) >= (data.summary.maxcoh.rev(i)*cutoff)));
end

% use the Steve method for cross-predictions

