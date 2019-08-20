%
%
%
%           by: david schoppik
%         date: 6/1/2006
%      purpose: to analyze the "trial-by-trial variability"
%               i.e. to generate and test the a linear filter
%               on the covariance between trials and spikes
%
%        usage: called by "prepare.m," this will return a whole set of
%               fields in the structure "unit" that it is passed:
%               unit = tbtv_final(unit,use_fraction) where use_fraction is
%               the % of trials to use for filter generation (0.4 is the
%               default, shown to minimize variance)
%
%               output fields include:
%               unit.filterbank <- the filters
%               unit.test <- the main analysis
%               unit.cross <- how does it do on the analysis of the mean
%               unit.bytrial <- the same analysis, broken down by trial


function [unit] = tbtv_final(unit,use_fraction)

if nargin == 1
  use_fraction = 0.4; % what percentage to use to generate the filter
end

binsize = 1;
numdraws = 150; % bootstrap parameter; how many draws to estimate the random filter
makefig = 1; % make a figure for each unit

% clear out the old values and preallocate
numtrials = size(unit.evel,1);
tlen = size(unit.binaryspikes,2);

unit.filterbank.var = single(zeros(tlen,numdraws));
unit.filterbank.randvar = single(zeros(tlen,numdraws));
unit.filterbank.best = single(zeros(tlen,numdraws));

unit.test.n = single(zeros(1,numdraws));
unit.test.corr = single(zeros(1,numdraws));
unit.test.timecorr = single(zeros(numdraws,tlen));
unit.test.mse = single(zeros(1,numdraws));
unit.test.timemse = single(zeros(numdraws,tlen));
unit.test.shufcorr = single(zeros(1,numdraws));
unit.test.shuftimecorr = single(zeros(numdraws,tlen));
unit.test.shufmsedist = single(zeros(1,numdraws));
unit.test.shuftimemse = single(zeros(numdraws,tlen));

unit.cross.corr = single(zeros(1,numdraws));
unit.cross.timecorr = single(zeros(numdraws,tlen));
unit.cross.bestcorr = single(zeros(1,numdraws));
unit.cross.shufcorr = single(zeros(1,numdraws));

unit.bytrial = nan(numdraws,numtrials);


%%%%%%%%%%%%%%%%%%%%%%
% general trial info %
%%%%%%%%%%%%%%%%%%%%%%

spikes = unit.binaryspikes;
eye = unit.evel;

% how many trials were there?

unit.general.numtrials = numtrials;

% some data on the firing rates
unit.general.fr.pref = mean(sum(spikes(:,201:end),2));
unit.general.fr.base = mean(sum(spikes(:,1:200),2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter generation and testing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% since we're bootstrapping, we first fill arrays first with randomly
% shuffled trials, so that we can use vector operations later
fulle3D = single(zeros([fliplr(size(eye)) numdraws]));
fullspikes3D = single(zeros([fliplr(size(spikes)) numdraws]));

% permute
[jnk perm] = sort(rand(numdraws,numtrials),2);

% fill the arrays
% I couldn't find an easy non-loop way to do this, though it no
% doubt exists with some tricky reshaping.
for draw = 1:numdraws
  fulle3D(:,:,draw) = eye(perm(draw,:),:)';
  fullspikes3D(:,:,draw) = spikes(perm(draw,:),:)';
end

% subtract the mean in time -- we're doing the analysis on the covariance
vare3D = fulle3D - repmat(mean(fulle3D,2),[1,numtrials,1]);
varspikes3D = fullspikes3D - repmat(mean(fullspikes3D,2),[1,numtrials,1]);

% split into filter and test trials
% find the index that corresponds to the use_fraction
stop = floor(use_fraction*numtrials);
useindex = (1:stop);
testindex = (stop+1:numtrials);
numfilttrials = length(useindex);
numtesttrials = length(testindex);

%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the matrices %
%%%%%%%%%%%%%%%%%%%%%%%%

% matrices of variability for filter generation
varE_filt = vare3D(:,useindex,:);
varS_filt = varspikes3D(:,useindex,:);

% matrices of variability for testing
varE_test = vare3D(:,testindex,:);
varS_test = varspikes3D(:,testindex,:);

% complete matrices for filter generation (subtract global mean)
fullE_filt = fulle3D(:,useindex,:) - mean(mean(mean(fulle3D(:,useindex,:))));
fullS_filt = fullspikes3D(:,useindex,:) - mean(mean(mean(fullspikes3D(:,useindex,:))));

% complete matrices for testing (subtract global mean)
fullE_test = fulle3D(:,testindex,:) - mean(mean(mean(fulle3D(:,testindex,:))));
fullS_test = fullspikes3D(:,testindex,:) - mean(mean(mean(fullspikes3D(:,testindex,:))));

% randomly permute the matrix of variability for filter generation to
% test the statistical significance of the filter
[jnk filtperm] = sort(rand(numdraws,numfilttrials),2);
for draw = 1:numdraws
  varS_rand(:,:,draw) = varS_filt(:,filtperm(draw,:),draw);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike-Triggered Averaging %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate the filters
varsta = getfilter(varS_filt,varE_filt);
shufvarsta = getfilter(varS_rand,varE_filt);
beststa = getfilter(fullS_filt,fullE_filt);

% save the filters
unit.filterbank.var = single(squeeze(varsta));
unit.filterbank.randvar = single(squeeze(shufvarsta));
unit.filterbank.best = single(squeeze(beststa));

%%%%%%%%%%%%%%%%%%%%%%
% Test the var model %
%%%%%%%%%%%%%%%%%%%%%%

% use the filters to make predictions about the variability on each trial
varE_pred = myconv(varS_test,repmat(varsta,[1,size(varS_test,2),1]));

fullE_pred = myconv(squeeze(mean(fullS_test,2)),squeeze(varsta));
bestE_pred = myconv(squeeze(mean(fullS_test,2)),squeeze(beststa));
shufE_pred = myconv(squeeze(mean(fullS_test,2)),squeeze(shufvarsta));



% test each filter against the shuffled data
for draw = 1:numdraws

  predicted_var = varE_pred(:,:,draw); % predicted
  actual_var = varE_test(:,:,draw); % actual

  % get some basic parameters on the relationship
  unit.test.n(draw) = single(length(predicted_var(:)));

  % analyze the relationship
  % (correlation coefficient and mean square error alone and in time)
  unit.test.corr(draw) = single(corr2(predicted_var(:),actual_var(:)));
  unit.test.mse(draw) = single(mean(mean((actual_var - predicted_var).^2,2)));
  unit.test.timecorr(draw,:) = single(mycorr(predicted_var',actual_var'));
  unit.test.timemse(draw,:) = single(mean((actual_var - predicted_var).^2,2));

  % measure how well each filter did at predicting the mean
  actual_full = single(mean(fullE_test(:,:,draw),2));
  predicted_full = single(fullE_pred(:,draw));
  shuf_predicted_full = single(shufE_pred(:,draw));
  best_predicted_full = single(bestE_pred(:,draw));

  unit.cross.corr(draw) = single(corr2(predicted_full(:),actual_full(:)));
  unit.cross.bestcorr(draw) = single(corr2(best_predicted_full(:),actual_full(:)));
  unit.cross.shufcorr(draw) = single(corr2(shuf_predicted_full(:),actual_full(:)));
  unit.cross.timecorr(draw,:) = single(mycorr(predicted_var',actual_var'));

  % bootstrap distributions for each analysis, for statistical comparison
  %for shuffledraw = 1:numshuffledraws
  shuffle_predicted_var = single(varE_pred(:,randperm(numtesttrials)));
  unit.test.shufcorr(draw) = single(corr2(shuffle_predicted_var(:),actual_var(:)));
  unit.test.shufmsedist(draw) = single(mean(mean((actual_var - shuffle_predicted_var).^2,2)));
  unit.test.shuftimecorr(draw,:) = single(mycorr(shuffle_predicted_var',actual_var'));
  unit.test.shuftimemse(draw,:) = single(mean((actual_var - shuffle_predicted_var).^2,2));
  unit.test.shufmse(draw) = single(mean(mean((actual_var - shuffle_predicted_var).^2,2)));
  %end

  % keep track of how we do for each trial
  unit.bytrial(draw,perm(draw,testindex)) = mycorr(predicted_var,actual_var);
  
end


%%%%%%%%%%%%%%%%%
% make a figure %
%%%%%%%%%%%%%%%%%


if makefig
  % get the sizing correct
  triallen = (size(eye,2) - 1) - 200;
  filtlen = (size(eye,2) - 1)/2;
  
  
  % make a figure
  h = figure;
  set(gcf,'position',[800 60 450 675],'paperpositionmode','auto','Inverthardcopy','on')

  % plot eye movement
  subplot(3,2,1)
  plot(-200:triallen,eye','k')
  set(gca,'xlim',[-200 triallen],'xtick',[-200 0 250 500 750]);
  title(['Unit ID: ',num2str(unit.unitid)])
  set(gca,'tickdir','out',...
    'box','off')
  xlabel('Time (ms)')
  ylabel('Eye velocity (deg/s)')

  % now plot the rasters
  subplot(3,2,2)
  hold on
  for i = 1:size(spikes,1)
    plot([find(spikes(i,:));find(spikes(i,:))],[i*ones(1,length(find(spikes(i,:))));(i+1)*ones(1,length(find(spikes(i,:))))],'k');
  end
  set(gca,'xlim',[1 size(eye,2)],'xtick',[1 200 450 700 950],'xticklabel',[-200 0 250 500 750]);
  set(gca,'ylim',[1 size(spikes,1)+1]);
 % title([unit.dir(end-7:end),'.',num2str(unit.localid)]);
  set(gca,'tickdir','out',...
    'box','off')
  xlabel('Time (ms)')

  % plot the filter
  subplot(3,2,3)
  hold on
  plot(-filtlen:filtlen,flipud(squeeze(unit.filterbank.randvar(:,:))),'color',[.6 .6 .6])
  plot(-filtlen:filtlen,flipud(squeeze(unit.filterbank.var(:,:))),'color',pretty('i'))
  set(gca,'tickdir','out',...
    'box','off','xlim',[-filtlen filtlen],'xtick',[-filtlen 0 filtlen])
  xlabel('Lag (ms)')
  ylabel('deg/sec/Spike')

  % plot overall performance
  corrvals = unit.test.corr(:);
  corrdist = unit.test.shufcorr(:,:);
  corrdist = corrdist(:);

  subplot(3,2,4)
  hold on
  [b h] = hist(corrdist);
  bar(h,b./sum(b),.9)
  patch = findobj(gca,'Type','patch');
  set(patch(1),'facecolor',[.6 .6 .6],'edgecolor','k')

  [b h] = hist(corrvals);
  bar(h,b./sum(b),.9)
  patch = findobj(gca,'Type','patch');
  set(patch(1),'facecolor',pretty('i'),'edgecolor','k')
  title('Variability Data')
  set(gca,'tickdir','out',...
    'box','off')
  set(gca,'xlim',[-.75 .75])
  xlabel('R')

  % plot overall performance
  crossvals = unit.cross.corr(:);
  bestcross = unit.cross.bestcorr(:);
  shufcross = unit.cross.shufcorr(:);

  subplot(3,2,5)
  set(gca,'visible','off')
  text(0,.9,sprintf('Cross: %0.1g (%0.1g)',mean(crossvals),std(crossvals)))
  text(0,.8,sprintf('Best Cross: %0.1g (%0.1g)',mean(bestcross),std(bestcross)))
  text(0,.7,sprintf('Shuffle Cross: %0.1g (%0.1g)',mean(shufcross),std(shufcross)))

  % plot the time course
  subplot(3,2,6)
  hold on
  plot(-200:triallen,squeeze(unit.test.timecorr(:,:)),'color',pretty('i'))
  shufvals = squeeze(unit.test.shuftimecorr(:,:,:));
  shufmean = squeeze(mean(shufvals));
  shufstd = 2.*squeeze(std(shufvals,0,1));
  plot(-200:triallen,shufmean,'color',[.6 .6 .6]);
  plot(-200:triallen,shufmean + shufstd,'--','color',[.6 .6 .6],'linewidth',3);
  plot(-200:triallen,shufmean - shufstd,'--','color',[.6 .6 .6],'linewidth',3);
  set(gca,'tickdir','out',...
    'box','off')
  set(gca,'xlim',[-200 triallen],'xtick',[-200 0 250 500 750]);
  title('Correlation in time')
  xlabel('Time (ms)')
  ylabel('R')


  %print('-dpsc2',h)
  %close(h)

end

end % main function

% use the real function, not the one here
% function r = getfilter(A,B)
% % A is the spikes
% % B is the eye movement
% 
% % prepare for the Fourier transforms
% nfft = 2^(1+nextpow2(size(A,1)));
% % to get rid of the edges from the xcorr/fft
% ftind = nfft/2-floor(size(A,1)/2)+1:nfft/2+round(size(A,1)/2);
% 
% % cross correlate
% % the conj is around A to flip the filter correctly
% STA = real(ifft(fft(A,nfft).*conj(fft(B,nfft))));
% 
% % get the power spectrum
% PS = mean(fft(A,nfft).*conj(fft(A,nfft)),2);
% 
% % decorrelate by division
% dc = fftshift(real(ifft(fft(mean(STA,2),nfft)./PS)),1);
% 
% % return the piece we want
% r = dc(ftind,:,:);
% 
% end

function r = myconv(A,B);

% prepare for the Fourier transforms
nfft = 2^(1+nextpow2(size(A,1)));

% to get rid of the edges from the xcorr/fft
ftind = nfft/2-floor(size(A,1)/2)+1:nfft/2+round(size(A,1)/2);

% transform and reshape
r = fftshift(real(ifft(fft(A,nfft).*conj(fft(B,nfft)))),1);

% resize
r = r(ftind,:,:);

end

function r = mycorr(A,B);
% not expecting things to be flipped
% i.e. time is across rows
% will only return the correlation along the diagonal
% which is identical to the unit of corr;

% get rid of the means
A = A-repmat(mean(A),size(A,1),1);
B = B-repmat(mean(B),size(B,1),1);

% compute the covariance
num = diag((A'*B));
denom = sqrt(diag(A'*A).*diag(B'*B));

r = num./denom;

end

%     hold on
%     [b h] = hist(shufcross);
%     bar(h,b./sum(b),.9)
%     patch = findobj(gca,'Type','patch');
%     set(patch(1),'facecolor',[.6 .6 .6],'edgecolor','k')
%
%     [b h] = hist(crossvals);
%     bar(h,b./sum(b),.9)
%     patch = findobj(gca,'Type','patch');
%     set(patch(1),'facecolor',pretty('i'),'edgecolor','k')
%
%     [b h] = hist(bestcross);
%     bar(h,b./sum(b),.9)
%     patch = findobj(gca,'Type','patch');
%     set(patch(1),'facecolor',pretty('a'),'edgecolor','k')
%
%     title('Variability Data')
%     set(gca,'tickdir','out',...
%       'box','off')
%     set(gca,'xlim',[-1 1])
%     xlabel('R')



%     subplot(3,2,6)
%     hold on
%     hist(unit.cross.randcorr(ind,:));
%     patch = findobj(gca,'Type','patch');
%     set(patch(1),'facecolor',[.6 .6 .6],'edgecolor','k')
%
%     hist(unit.cross.corr(ind,:));
%     patch = findobj(gca,'Type','patch');
%     set(patch(1),'facecolor',pretty('t'),'edgecolor','k')
%     title('All data')
%
%     set(gca,'tickdir','out',...
%       'box','off','xlim',[-1 1])
%     xlabel('R')
%

%     % what's the dot product of the normalized filters (a measure of
%     % similarity of the shapes?)
%     unit.cross.dp(ind,draw) = (varsta/norm(varsta))*(fullsta/norm(fullsta))';
%
%     full_filter_on_variance = myconv(testspikes_in_time',repmat(fullsta,size(testspikes_in_time,1),1)');
%     var_filter_on_full = myconv(testspikes',repmat(varsta,size(testspikes,1),1)');
%
%     f_on_v = full_filter_on_variance(:);
%     v_on_f = var_filter_on_full(:);
%
%     unit.cross.var.corr.fonv(ind,draw) = corr2(f_on_v,vary);
%     unit.cross.var.corr.vonf(ind,draw) = corr2(v_on_f,y);
%
%   end



%   % we can only test the "variance" filter against the variance in the
%   % eye movement (i.e. subtract the mean in time
%   teste_in_time = fulle(testindex,:) - repmat(mean(fulle(testindex,:)),size(fulle(testindex,:),1),1);
%   testspikes_in_time = fullspikes(testindex,:) - repmat(mean(fullspikes(testindex,:)),size(fullspikes(testindex,:),1),1);
%   testrandspikes_in_time = randfullspikes(testindex,:) - repmat(mean(randfullspikes(testindex,:)),size(randfullspikes(testindex,:),1),1);
%


%
%   randvarx = randvarehat(:);
%   unit.test.var.randscalefac(ind,draw,:) = polyfit(randvarx(useind),vary(useind),1);
%   unit.test.var.randcorr(ind,draw) = corr2(randvarx(useind),vary(useind));
%   unit.test.var.randmse(ind,draw) = mean((vary(useind)-polyval(unit.test.var.randscalefac(ind,draw,:),randvarx(useind))).^2);
%   unit.test.var.randtimecorr(ind,draw,:) = mycorr(randvarehat,teste_in_time);
%

%  %%%%%%%%%%%%%%%%%%%%%%%
%   % Test the full model %
%   %%%%%%%%%%%%%%%%%%%%%%%
%
%     % convolve the spikes and filter to generate a prediction of the eye
%     % movements
%     fullehat = myconv(testspikes',repmat(fullsta,size(testspikes,1),1)');
%     randehat = myconv(testrandspikes',repmat(randsta,size(testrandspikes,1),1)');
%
%     % bin the data to test the predictions
%     x = fullehat(:);
%     y = teste(:)+mean(mean(fulle(testindex,:))); % add the mean back in here
%
%     useind = (1:length(x));
%
%     % fit (and check the fit) to the line for the linear part of the I/O fn
%     % scale the prediction, and check the MSE
%     unit.test.full.scalefac(ind,draw,:) = polyfit(x(useind),y(useind),1);
%     unit.test.full.corr(ind,draw) = corr2(x(useind),y(useind));
%     unit.test.full.mse(ind,draw) = mean((y(useind)-polyval(unit.test.full.scalefac(ind,draw,:),x(useind))).^2);
%     unit.test.full.n(ind,draw) = length(useind);
%     unit.test.full.timecorr(ind,draw,:) = mycorr(fullehat,teste);
%
%     randx = randehat(:);
%     unit.test.full.randscalefac(ind,draw,:) = polyfit(randx(useind),y(useind),1);
%     unit.test.full.randcorr(ind,draw) = corr2(randx(useind),y(useind));
%     unit.test.full.randmse(ind,draw) = mean((y(useind)-polyval(unit.test.full.randscalefac(ind,draw,:),randx(useind))).^2);
%     unit.test.full.randtimecorr(ind,draw,:) = mycorr(randehat,teste);
%
%
%
%     % Let's only test against the remaining fraction of the trials
%     teste = fulle(testindex,:) - mean(mean(fulle(testindex,:)));
%     testspikes = fullspikes(testindex,:) - mean(mean(fullspikes(testindex,:)));
%     testrandspikes = randfullspikes(testindex,:);


% function r = getfilternodec(A,B)
% % A is the spikes
% % B is the eye movement
%
%  % prepare for the Fourier transforms
% nfft = 2^(1+nextpow2(size(A,1)));
% % to get rid of the edges from the xcorr/fft
% ftind = nfft/2-floor(size(A,1)/2)+1:nfft/2+round(size(A,1)/2);
%
% % cross correlate
% STA = fftshift(real(ifft(fft(A,nfft).*conj(fft(B,nfft)))),1)';
%
% % return the piece we want
% r = mean(STA(:,ftind));
%
% end

%
%
%
%
%
%
%   % first, remove the OVERALL mean values
%   evar = e-mean(mean(e));
%   spvar = spikes-mean(mean(spikes));
%   randspvar = randspikes-mean(mean(randspikes));
%
%      % and do it Bill's way to get the variability filter
%     evar_in_time = e-repmat(mean(e),size(e,1),1);
%     spvar_in_time = spikes-repmat(mean(spikes),size(spikes,1),1);
%     randspvar_in_time = randspikes-repmat(mean(randspikes),size(spikes,1),1);
%
%     fullsta = getfilter(spvar',evar');
%     varsta = getfilter(spvar_in_time',evar_in_time');
%
%     randsta = getfilter(randspvar',evar');
%     randvarsta = getfilter(randspvar_in_time',evar_in_time');
%
%
%
%     % save the filters
%     unit.filterbank.var(ind,:,draw) = varsta;
%     unit.filterbank.randvar(ind,:,draw) = randvarsta;
%     unit.filterbank.full(ind,:,draw) = fullsta;
%     unit.filterbank.randfull(ind,:,draw) = randsta;

%
%
%    for draw = 1:numdraws
%     % randomly select a percentage of some of the trials
%     index = randperm(size(fulle,1));
%     stop = floor(use_fraction*size(fulle,1));
%     useindex = index(1:stop);
%     testindex= index(stop+1:end);
%
%     e = fulle(useindex,:);
%     spikes = fullspikes(useindex,:);
%
%     % generate an index vector to shuffle the trials that we are using
%     shufeye = randperm(size(e,1));
%     shuftestindex = testindex(randperm(length(testindex)));
%
%     % randomly distribute the spikes
%     %randfullspikes = reshape(fullspikes(randperm(prod(size(fullspikes)))),size(fullspikes,1),size(fullspikes,2));
%     randfullspikes = fullspikes(randperm(size(fullspikes,1)),:);
%     randspikes = randfullspikes(useindex,:);
%     randfr = reshape(fr(randperm(prod(size(fr)))),size(fr,1),size(fr,2));
%
%     % prepare for the Fourier transforms
%     nfft = 2^nextpow2(2*size(fulle,2)-1);
%
%     % to get rid of the edges from the xcorr/fft
%     goodind = (size(e,2)/2):(2*size(e,2)-1)-(size(e,2)/2);
%     ftind = nfft/2-floor(size(e,2)/2)+1:nfft/2+round(size(e,2)/2);


% only look at the variance in the bins that reflect from a percentage of
% the values for velocity found in the trials
% high = .85*max(mean(fulle(testindex,:)));
% low = .15*max(mean(fulle(testindex,:)));
% useind = find(y > low & y < high);

% now, remove the part of the correlation that is not due to the
% trials:
% meansta = getfilter(spvar',repmat(mean(evar),size(evar,1),1)');
% varsta = fullsta-meansta;
% randvarsta = randsta-randmeansta;
% randmeansta = getfilter(randspvar',repmat(mean(evar),size(evar,1),1)');
% unit.filterbank.mean(ind,:,draw) = meansta;
% unit.filterbank.bill(ind,:,draw) = billter;

%
% function r = mypowerspectrum(A)
% % takes an array of spikes and returns a vector corresponding to the power
% % spectrum -- returns in Fourier space appropriate for the denominator
% nfft = 2^(1+nextpow2(size(A,1)));
%
% F = fft(A,nfft);
% F = F.*conj(F);
%
% r = mean(F,2);
%
% end
%
% function r = mydecorrelate(a,b);
% % a should be a filter
% % b should be a power spectrum
%
% % prepare for the Fourier transforms
% nfft = 2^(1+nextpow2(size(a,1)));
% ftind = nfft/2-floor(size(a,1)/2)+1:nfft/2+round(size(a,1)/2);
%
% r = (real(ifft(fft(a,nfft)./b)));
%
% % resize
% r = r(ftind)';
%
% end
%
%
%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the heatmap and its mask,make a population version %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%warning off
%[r p] = corr(fr,eye);
%warning on
%mask = zeros(size(p));
%mask(find(p < .05/prod(size(p)))) = 1;
%unit.general.popheatmap = unit.general.popheatmap + mask;

% get the best point
%[frind eind] = find(abs(r) == max(max(abs(r))));



% % now smooth the filter by the number of neurons that actually went
% % into each point.  there should be a better way to do this....
% [jnk spikeindices] = find(spikes);
% smoothingfilt = zeros(1,2*size(spikes,2)-1);
% smoothingfilt(size(spikes,2)) = length(find(spikes));
% preinds = 1 - spikeindices;
% postinds = size(spikes,2) - spikeindices;
% for i = 1:round(size(spikes,2)/4)-1
%   smoothingfilt(i+round(size(spikes,2)/4)) = length(find(postinds >= i));
%   smoothingfilt(-i+round(size(spikes,2)/4)) = length(find(preinds <= -i));
% end
% smoothingfilt = smoothingfilt./length(find(spikes));
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% process the many draws %
%%%%%%%%%%%%%%%%%%%%%%%%%%


%unit.test.sta.best.r(ind,draw) = corr2(epred(:),teste(:));
%unit.test.sta.rand.r(ind,draw) = corr2(randepred(:),teste(:));
%unit.test.sta.shuf.r(ind,draw) = corr2(shufepred(:),shufteste(:));


% calculate the means and std for the resampled data
%   unit.test.best.mean(ind) = nanmean(unit.test.best.goffit(ind,:));
%   unit.test.best.std(ind) = nanstd(unit.test.best.goffit(ind,:));
%   unit.test.rand.mean(ind) = nanmean(unit.test.rand.goffit(ind,:));
%   unit.test.rand.std(ind) = nanstd(unit.test.rand.goffit(ind,:));
%   unit.test.shuf.mean(ind) = nanmean(unit.test.shuf.goffit(ind,:));
%   unit.test.shuf.std(ind) = nanstd(unit.test.shuf.goffit(ind,:));

%   unit.sig.rand(ind) = ttest(unit.test.best.goffit(ind,:),unit.test.rand.goffit(ind,:),.05/numdraws);
%   unit.sig.shuf(ind) = ttest(unit.test.best.goffit(ind,:),unit.test.shuf.goffit(ind,:),.05/numdraws);


%%%%%%%%%%
% old IO %
%%%%%%%%%%

%     bins = min(x):(max(x)-min(x))/25:max(x);
%
%     meanvals = [];
%     stdvals = [];
%     for i = 1:length(bins)-1;
%       stdvals(i) = nanstd(y(find(x > bins(i) & x < bins(i+1))));
%       meanvals(i) = nanmean(y(find(x > bins(i) & x < bins(i+1))));
%     end
%
%
%
%     if ~isempty(roi)
%       % we can fit a line here to scale the filter
%       unit.io.scalefac{ind,draw} = polyfit(bins(roi),meanvals(roi),1);
%       % check to see that it is a good fit
%       unit.io.corr{ind,draw} = corr2(meanvals(roi),polyval(unit.io.scalefac{ind,draw},bins(roi)));
%       % keep track of the (r)egion (o)f (i)nterest
%       unit.io.roi{ind,draw} = roi;
%
%       % and only calculate the variance from the roi
%       drawnstd(draw) = nanmean(stdvals(roi));
%     else
%       % we can fit a line here to scale the filter
%       unit.io.scalefac{ind,draw} = nan;
%       % check to see that it is a good fit
%       unit.io.corr{ind,draw} = nan;
%       % keep track of the (r)egion (o)f (i)nterest
%       unit.io.roi{ind,draw} = nan;
%
%       % and only calculate the variance from the roi
%       drawnstd(draw) = nan;
%     end
%
%     unit.io.meanvals{ind,draw} = meanvals;
%     unit.io.stdvals{ind,draw} = stdvals;
%     unit.io.bins{ind,draw} = bins(1:end-1);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Moment-by-moment analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       % get the heatmap and its mask,make a population version
%       % since we will estimate the eye movement at the end by convolving with
%       % the firing rate, put firing rate first here so that the filters come out facing
%       % the right way.
%       warning off
%       [mbm] = corr(eye(useindex,:),fr(useindex,:));
%       [randmbm] = corr(eye(useindex,:),randfr(useindex,:));
%       [shufmbm] = corr(eye(shufeye,:),fr(useindex,:));
%       warning on
%
%       % get the slope at the best point (determined once per experiment)
%       [beta] = polyfit(fr(useindex,frind(1)),eye(useindex,eind(1)),1);
%       unit.mbm.best.slope(ind,draw) = beta(1);
%
%       [shufbeta] = polyfit(fr(useindex,frind(1)),eye(shufeye,eind(1)),1);
%       unit.mbm.shuf.slope(ind,draw) = shufbeta(1);
%
%       % get the MSE and R^2 on the model fit
%       % tested against the unused trials
%       unit.mbm.best.mse(ind,draw) = mean((eye(testindex,eind(1)) - polyval(beta,fr(testindex,frind(1)))).^2);
%       unit.mbm.best.goffit(ind,draw) = (corr2(eye(testindex,eind(1)),polyval(beta,fr(testindex,frind(1))))).^2;
%
%       unit.mbm.shuf.mse(ind,draw) = mean((eye(shuftestindex,eind(1)) - polyval(shufbeta,fr(testindex,frind(1)))).^2);
%       unit.mbm.shuf.goffit(ind,draw) = (corr2(eye(shuftestindex,eind(1)),polyval(shufbeta,fr(testindex,frind(1))))).^2;
%
%       % get the MBM "filter" via diagonalization
%       for i = -(length(mbm)-1):length(mbm)-1
%         mbmF1(i+length(mbm)) = nanmean(diag(mbm,i));
%         mbmrandF1(i+length(mbm)) = nanmean(diag(randmbm,i));
%         mbmshufF1(i+length(mbm)) = nanmean(diag(shufmbm,i));
%       end
%
%       mbmF1 = mbmF1(goodind);
%       mbmrandF1 = mbmrandF1(goodind);
%       mbmshufF1 = mbmshufF1(goodind);
%
%       unit.mbm.filterbank.best(ind,:,draw) = mbmF1;
%       unit.mbm.filterbank.rand(ind,:,draw) = mbmrandF1;
%       unit.mbm.filterbank.shuf(ind,:,draw) = mbmshufF1;
%
%       % test the prediction
%       % we're only testing against the remaining 20% of the trials
%       testeye = eye(testindex,:);
%       shuftesteye = eye(shuftestindex,:);
%       testfr = fr(testindex,:);
%       testrandfr = randfr(testindex,:);
%
%       eyehat = fftshift(real(ifft(fft(testfr',nfft).*conj(fft(repmat(mbmF1,size(testfr,1),1)',nfft)))))';
%       randeyehat = fftshift(real(ifft(fft(testfr',nfft).*conj(fft(repmat(mbmrandF1,size(testfr,1),1)',nfft)))))';
%       shufeyehat = fftshift(real(ifft(fft(testfr',nfft).*conj(fft(repmat(mbmshufF1,size(testfr,1),1)',nfft)))))';
%
%       eyepred = eyehat(:,ftind);
%       randeyepred = randeyehat(:,ftind);
%       shufeyepred = shufeyehat(:,ftind);
%
%       unit.mbm.filter.best.goffit(ind,draw) = corr2(eyepred(:),testeye(:)).^2;
%       unit.mbm.filter.shuf.goffit(ind,draw) = corr2(randeyepred(:),testeye(:)).^2;
%       unit.mbm.filter.shuf.goffit(ind,draw) = corr2(shufeyepred(:),shuftesteye(:)).^2;
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Spike-Triggered Cov?????? %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % we'll do a spike-triggered and a random-triggered average
%
%     % remove the time-based mean values
%     evar = e-repmat(mean(e),size(e,1),1);
%     spvar = spikes-repmat(mean(spikes),size(spikes,1),1);
%     randspvar = randspikes-mean(mean(randspikes));
%
%     xc = fftshift(real(ifft(fft(spvar',nfft).*conj(fft(evar',nfft)))))';
%     randxc = fftshift(real(ifft(fft(randspvar',nfft).*conj(fft(evar',nfft)))))';
%     shufxc = fftshift(real(ifft(fft(spvar',nfft).*conj(fft(evar(shufeye,:)',nfft)))))';
%
%     F1 = mean(xc(:,ftind))/norm(mean(xc(:,ftind)));
%     randF1 = mean(randxc(:,ftind))/norm(mean(randxc(:,ftind)));
%     shufF1 = mean(shufxc(:,ftind))/norm(mean(shufxc(:,ftind)));
%
%     % save the filters
%     unit.filterbank.stc.best(ind,:,draw) = F1;
%     unit.filterbank.stc.rand(ind,:,draw) = randF1;
%     unit.filterbank.stc.shuf(ind,:,draw) = shufF1;
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Test the filter predicitons %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     % we're only testing against the remaining fraction of the trials
%     teste = fulle(testindex,:);
%
%     shufteste = fulle(shuftestindex,:);
%
%     testspikes = fullspikes(testindex,:);
%     testrandspikes = randfullspikes(testindex,:);
%
%     ehat = fftshift(real(ifft(fft(testspikes',nfft).*conj(fft(repmat(F1,size(testspikes,1),1)',nfft)))))';
%     randehat = fftshift(real(ifft(fft(testrandspikes',nfft).*conj(fft(repmat(randF1,size(testrandspikes,1),1)',nfft)))))';
%     shufehat = fftshift(real(ifft(fft(testspikes',nfft).*conj(fft(repmat(shufF1,size(testspikes,1),1)',nfft)))))';
%
%     epred = ehat(:,ftind);
%     randepred = randehat(:,ftind);
%     shufepred = shufehat(:,ftind);
%
%     unit.test.stc.best.goffit(ind,draw) = corr2(epred(:),teste(:)).^2;
%     unit.test.stc.rand.goffit(ind,draw) = corr2(randepred(:),teste(:)).^2;
%     unit.test.stc.shuf.goffit(ind,draw) = corr2(shufepred(:),shufteste(:)).^2;
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Look at the I/O function %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     % first generate a prediction
%     ehat = fftshift(real(ifft(fft(spikes',nfft).*conj(fft(repmat(F1,size(spikes,1),1)',nfft)))))';
%     epred = ehat(:,ftind);
%
%     % bin the data to test the predictions
%     x = epred(:);
%     y = evar(:);%+(mean(mean(e)));
%     bins = min(x):(max(x)-min(x))/25:max(x);
%
%     meanvals = [];
%     stdvals = [];
%     for i = 1:length(bins)-1;
%       stdvals(i) = nanstd(y(find(x > bins(i) & x < bins(i+1))));
%       meanvals(i) = nanmean(y(find(x > bins(i) & x < bins(i+1))));
%     end
%
%     % only look at the variance in the bins that reflect from 10%-90% of
%     % the values for velocity found in the trials
%     high = .9*max(mean(fulle));
%     low = .1*max(mean(fulle));
%     highind = find(meanvals(1:end-1) < high & meanvals(2:end) > high);
%     lowind = find(meanvals(1:end-1) < low & meanvals(2:end) < low);
%     roi = (lowind:highind); % the linear part, we hope
%
%     if ~isempty(roi)
%       % we can fit a line here to scale the filter
%       unit.io.stc.scalefac{ind,:,draw} = polyfit(bins(roi),meanvals(roi),1);
%       % check to see that it is a good fit
%       unit.io.stc.corr{ind,:,draw} = corr2(meanvals(roi),polyval(unit.io.scalefac{ind,:,draw},bins(roi)));
%       % keep track of the (r)egion (o)f (i)nterest
%       unit.io.stc.roi{ind,:,draw} = roi;
%
%       % and only calculate the variance from the roi
%       drawnstd(draw) = nanmean(stdvals(roi));
%     else
%       % we can fit a line here to scale the filter
%       unit.io.stc.scalefac{ind,:,draw} = nan;
%       % check to see that it is a good fit
%       unit.io.stc.corr{ind,:,draw} = nan;
%       % keep track of the (r)egion (o)f (i)nterest
%       unit.io.stc.roi{ind,:,draw} = nan;
%
%       % and only calculate the variance from the roi
%       drawnstd(draw) = nan;
%     end
%
%     unit.io.stc.meanvals{ind,:,draw} = meanvals;
%     unit.io.stc.stdvals{ind,:,draw} = stdvals;
%     unit.io.stc.bins{ind,:,draw} = bins(1:end-1);

