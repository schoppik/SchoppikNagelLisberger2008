% clean up from previous runs
close all
clear all

% a function to simulate the consequences of modeling the residual eye
% movements against a non-stationary mean
% makes a supplementary figure for tbtv

% variables
numtrials = 100;
maxFR = .50;
baselineFR = 0.005;
S1 = figure;
numdraws = 2;
numtrialstodraw = 100;

% make a filter
F = (gaussline([0 100 5],-350:350));
%F = [zeros(1,100) gaussian([0 35 1],-300:300)];
F = F./norm(F);

% get the functions to scale the residuals from actual data
cd /Users/schoppik/Documents/Data/
load example
realEye = example.evel;
realEyeResiduals = realEye - repmat(mean(realEye),[95 1]);

% the smoothing function is the square of the autocorrelation function of
% the real data
smfun = (mean(myconv(realEyeResiduals',realEyeResiduals'),2));
smfun = smfun.^2;
% increase the size of the integral of the filter so that it preserves the
% original values
%smfun = (smfun./sum(smfun)).*12.5;
smfun = smfun./norm(smfun);

% generate some residual eye movements
residualEye = randn(numtrials,701);

% smooth the residuals
residualEye = myconv(residualEye',(repmat((smfun)',[numtrials 1])'))';

% scale the residuals
residualEye = repmat(std(realEyeResiduals),[numtrials 1]).*residualEye;

% make an average eye velocity trace
meanEye = mean(realEye);

subplot(2,3,4)
hold on
plot(-200:500,residualEye','color',[.6 .6 .6])
plot(-200:500,std(residualEye),'k','linewidth',2)
set(gca,'tickdir','out','xlim',[-200 500],'ylim',[-10 10])
xlabel('Time (ms)')
ylabel('Residual velocity (deg/s)')

% and now the mean eye movements
fullEye = residualEye + repmat(meanEye,[numtrials 1]);
subplot(2,3,1)
plot(-200:500,fullEye','color',[.6 .6 .6])
set(gca,'tickdir','out','xlim',[-200 500],'box','off','ylim',[-10 30])
xlabel('Time (ms)')
ylabel('Velocity (deg/s)')

% now get the firing rate from the mean
firingRate = (myconv(fullEye',repmat(F,[numtrials 1])'))';

% rectify
firingRate(find(firingRate<0)) = 0;

% scale so we get a reasonable number of spikes and a baseline
firingRate = maxFR*(firingRate./max(max(firingRate)));
firingRate = baselineFR + firingRate;

% now translate that into spikes 
spikes = zeros(numtrials,701);
randomComparator = rand(numtrials,701);
spikes(find(randomComparator < firingRate)) = 1;

subplot(2,3,2)
[x,y] = quickraster_phy(spikes,200,500);
plot(x,y,'k')
set(gca,'xlim',[-200 500],'ylim',[0 100],'box','off')

% get the residual spikes
residualSpikes = spikes - repmat(mean(spikes),[numtrials 1]);
subplot(2,3,5)
hold on
plot(-200:500,residualSpikes','color',[.6 .6 .6])
plot(-200:500,mean(residualSpikes) + std(residualSpikes),'color','k','linewidth',2)
set(gca,'xlim',[-200 500],'ylim',[-0.25 1])

for i = 1:numdraws
  dex = randperm(numtrials);
  Fhat(i,:) = flipud(getfilter(residualSpikes(dex(1:numtrialstodraw),:)',...
    residualEye(dex(1:numtrialstodraw),:)'));
  fullFhat(i,:) = flipud(getfilter(spikes(dex(1:numtrialstodraw),:)',...
    fullEye(1:numtrialstodraw,:)'));

  Fhat(i,:) = Fhat(i,:)./norm(Fhat(i,:));
  fullFhat(i,:) = fullFhat(i,:)./norm(fullFhat(i,:)); 
end

% plot the filters
subplot(2,3,3)
hold on
errorbar(-350:350,mean(fullFhat),std(fullFhat),'color','k')
plot(-350:350,F./norm(F),'color',[.6 .6 .6],'linewidth',2)
set(gca,'xlim',[-350 350],'box','off','tickdir','out')
xlabel('Time (ms)')
ylabel('Normalized units')

subplot(2,3,6)
hold on
errorbar(-350:350,mean(Fhat),std(Fhat),'color','k');
plot(-350:350,F./norm(F),'color',[.6 .6 .6],'linewidth',2)
set(gca,'xlim',[-350 350])
xlabel('Time (ms)')
ylabel('Normalized units')
