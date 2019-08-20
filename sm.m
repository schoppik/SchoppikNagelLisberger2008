%%%%%%%%%%%%%
% Figure S1 %
%%%%%%%%%%%%%
% describing the alignment and the position data

fS1 = figure;
subplot(1,3,1)
hold on
plot([0 1],[0 1],'k','linewidth',2)
plot(data.summary.maxcorr.rev(data.sigind),data.summary.maxcorr.pos(data.sigind),'kx','markersize',6)
plot(data.summary.maxcorr.rev(data.nsigind),data.summary.maxcorr.pos(data.nsigind),'.',...
  'markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','xlim',[0 1],'xtick',[0 .5 1],'ylim',[0 1],'ytick',[0 .5 1])
title(sprintf('r = %0.2g',corr2(data.summary.maxcorr.rev(data.sigind),data.summary.maxcorr.pos(data.sigind))));
ylabel('Position Trace')

subplot(1,3,2)
hold on
plot([0 1],[0 1],'k','linewidth',2)
plot(data.summary.maxcorr.rev(data.sigind),data.summary.maxcorr.aln(data.sigind),'kx','markersize',6)
plot(data.summary.maxcorr.rev(data.nsigind),data.summary.maxcorr.aln(data.nsigind),'.',...
  'markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','xlim',[0 1],'xtick',[0 .5 1],'ylim',[0 1],'ytick',[0 .5 1])
title(sprintf('r = %0.2g',corr2(data.summary.maxcorr.rev(data.sigind),data.summary.maxcorr.aln(data.sigind))));
ylabel('Aligned Trace')
xlabel('Velocity Trace')

subplot(1,3,3)
hold on
plot([0 1],[0 1],'k','linewidth',2)
plot(data.summary.maxcorr.rev(data.sigind),data.summary.maxcorr.srt(data.sigind),'kx','markersize',6)
plot(data.summary.maxcorr.rev(data.nsigind),data.summary.maxcorr.srt(data.nsigind),'.',...
  'markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','xlim',[0 1],'xtick',[0 .5 1],'ylim',[0 1],'ytick',[0 .5 1])
title(sprintf('r = %0.2g',corr2(data.summary.maxcorr.rev(data.sigind),data.summary.maxcorr.srt(data.sigind))));
ylabel('Short Trace')

fie(fS1,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS1.xml')
print(fS1,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS1.eps')

%%%%%%%%%%%%%
% Figure S2 %
%%%%%%%%%%%%%
% figure describing the consequences of binning

dex = find(example.alltype == example.prefdir);
nfft = 2048;
f = [0:nfft-1]*1000/nfft;
windowing = [1 20];
evar = example.allevel(dex,:) - repmat(mean(example.allevel(dex,:)),[length(dex) 1]);


for i = 1:length(windowing)
  spvar(:,:,i) = windowspikedata(example.allspikes(dex,:),windowing(i)) - repmat(mean(windowspikedata(example.allspikes(dex,:),windowing(i))),[length(dex) 1]);
   % get the power spectrum
   pse(i,:) = mean(abs(fft(spvar(:,:,i)',nfft).*conj(fft(evar',nfft))),2)+eps;
   ps(i,:) = mean(fft(spvar(:,:,i)',nfft).*conj(fft(spvar(:,:,i)',nfft)),2)+eps;
  D(i,:) = getfilter(spvar(:,:,i)',evar');
end

fS2 = figure;
subplot(2,3,1)
plot(-200:500,spvar(:,:,1)','color',[.6 .6 .6])
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200 0 250 500],'ylim',[-0.25 1],'ytick',[-0.25 0 0.5 1],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Time (ms)')
ylabel('Spikes/bin')

subplot(2,3,2)
plot(-350:350,fliplr(D(1,:)),'k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-350 350],'xtick',[-350:350:350],'ylim',[-0.5 1],'ytick',[-0.5:.5:1],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Lag (ms)')
ylabel('Deg/s/spike')
title('1 ms bin size')

fudgefact = 0.005/0.04;
subplot(2,3,3)
hold on
bar(f(1:nfft/2),ps(1,1:(nfft/2))./sum(ps(1,(1:nfft/2))),'facecolor','k')
plot(f(1:nfft/2),pse(1,1:nfft/2)./sum(pse(1:nfft/2)).*fudgefact,'r','linewidth',3)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[0 100],'xtick',[0:25:100],'ylim',[0 .005],'ytick',[],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
ylabel('Power (norm)')
xlabel('Frequency (Hz)')

subplot(2,3,4)
plot(-200:500,spvar(:,:,2)','color',[.6 .6 .6])
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200 0 250 500],'ylim',[-0.25 0.25],'ytick',[-0.25:0.25:0.25],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Time (ms)')
ylabel('Spikes/bin')

subplot(2,3,5)
plot(-350:350,fliplr(D(2,:)),'k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-350 350],'xtick',[-350:350:350],'ylim',[-0.5 1],'ytick',[-0.5:.5:1],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Lag (ms)')
ylabel('Deg/s/spike')
title('20ms bin size')

fudgefact = 0.025/0.06;
subplot(2,3,6)
hold on
bar(f(1:nfft/2),ps(2,1:nfft/2)./sum(ps(2,1:nfft/2)),'facecolor','k')
plot(f(1:nfft/2),pse(2,1:nfft/2)./sum(pse(1:nfft/2)).*fudgefact,'r','linewidth',3)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[0 100],'xtick',[0:25:100],'ylim',[0 .025],'ytick',[],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)

xlabel('Frequency')
ylabel('Power (norm)')

fie(fS2,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS2.xml')
print(fS2,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS2.eps')

%%%%%%%%%%%%%
% Figure S3 %
%%%%%%%%%%%%%

fS3a = figure;

subplot(1,2,1)
hold on
plot(x(randdex(:)),y(randdex(:)),'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',3)
errorbar(bins,meanvals,stdvals,stdvals,'k.','markersize',18);
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-2.5 2.5],'xtick',[-2.5:2.5:2.5],'ylim',[-7.5 7.5],'ytick',[-7.5:7.5:7.5],'tickdir','out',...
   'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Predicted variability (deg/s)');
ylabel('Actual variability (deg/s)');
[r p] = corrcoef(x(:),y(:));
title(['r  = ',num2str(r(1,2)),' n = ',num2str(length(x)),' p = ',num2str(p(1,2))])

subplot(1,2,2)
hold on
hist(example.test.shufcorr,[-.5:.1:.5])
hist(example.test.corr,[.2:.05:.6])
h = findobj(gca,'type','patch');
set(h(1),'facecolor','k','edgecolor',[.6 .6 .6])
set(h(2),'facecolor',[.6 .6 .6],'edgecolor','k')
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'ylim',[0 80],'ytick',[0:20:80],...
  'xlim',[-.8 .8],'xtick',[-.8:.4:.8],...
   'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('R')
ylabel('Draws')


fie(fS3a,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS3a.xml')
print(fS3a,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS3a.eps')


% use the filter to generate some predictions
spvar = example.binaryspikes - repmat(mean(example.binaryspikes),[size(example.binaryspikes,1),1]);
evar = example.evel - repmat(mean(example.evel),[size(example.evel,1),1]);
ehat = myconv(spvar',repmat(mean(example.filterbank.var,2),[1,size(example.binaryspikes,1)]))';

m = nanmean(example.bytrial);
s = nanstd(example.bytrial);
[jnk dex] = sort(m);

worstid = dex(1);
bestid = dex(end);
midid = dex(round(length(dex)/2));

x = ehat(:);
y = evar(:);

randdex = randperm(length(x));

[jnk xind] = sort(x);

numbins = 25;

numperbin = floor(size(x,1)/numbins);

meanvals = [];
stdvals = [];

bins = [];
for i = 1:numbins
  meanvals(i) = mean(y(xind(((i-1)*numperbin)+1:i*numperbin)));
  stdvals(i) = std(y(xind(((i-1)*numperbin)+1:i*numperbin)));
  bins(i) = mean(x(xind(((i-1)*numperbin)+1:i*numperbin)));
end

tc = squeeze(example.test.timecorr);
rand_tc = squeeze(example.test.shuftimecorr);
mask = nan(1,701);
sig_time = ttest(tc,rand_tc,0.05/250);
mask(find(sig_time)) = 1;

fS3b = figure;
subplot(1,3,1)
hold on
plot(-200:500,2+evar(bestid,:)./max(abs(evar(bestid,:))),'k--','linewidth',2)
plot(-200:500,2+ehat(bestid,:)./max(abs(ehat(bestid,:))),'k','linewidth',2)
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'tickdir','out',...
   'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
ylabel('Eye variability (norm)');
xlabel('Time (ms)')
title(sprintf('r = %0.2g (%0.2g)',m(bestid),s(bestid)))  

subplot(1,3,2)
hold on
plot(-200:500,evar(worstid,:)./max(abs(evar(midid,:))),'k--','linewidth',2)
plot(-200:500,ehat(worstid,:)./max(abs(ehat(midid,:))),'k','linewidth',2)
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'tickdir','out',...
   'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
ylabel('Eye variability (norm)');
xlabel('Time (ms)')
title(sprintf('r = %0.2g (%0.2g)',m(midid),s(midid)))  

subplot(1,3,3)
hold on
errorbar(1:95,m(dex),s(dex),s(dex),'k')
hline(2.*(mean(example.test.shufcorr) + std(example.test.shufcorr)))
hline(2.*(mean(example.test.shufcorr) - std(example.test.shufcorr)))
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[0 100],'xtick',[0:25:100],'tickdir','out',...
   'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
ylabel('r')
xlabel('Trial number')

fie(fS3b,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS3b.xml')
print(fS3b,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS3b.eps')


%%%%%%%%%%%%%
% Figure S4 %
%%%%%%%%%%%%%

fS4 = figure;

subplot(3,2,1)
hold on
hist(data.crosscorr.rev(data.sigind),[-1:.2:1]);
hist(data.crosscorr.rev(data.nsigind),[-1:.2:1]);
h = findobj(gca,'type','patch');
set(h(2),'facecolor','k','edgecolor',[.6 .6 .6])
set(h(1),'facecolor',[.6 .6 .6],'edgecolor','k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-1.1 1.1],'xtick',[-1:0.5:1],'ylim',[0 25],'ytick',[0:12.5:25],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Value of r')
ylabel('Number of neurons')

%indbest = ((data.bestcrosscorr.rev - abs(data.crosscorr.rev)) ./ (data.bestcrosscorr.rev + abs(data.crosscorr.rev)));
indbest = data.crosscorr.rev./data.bestcrosscorr.rev;

subplot(3,2,2)
hold on
hist(indbest(data.sigind),[-1:.2:1]);
hist(indbest(data.nsigind),[-1:.2:1]);
h = findobj(gca,'type','patch');
set(h(2),'facecolor','k','edgecolor',[.6 .6 .6])
set(h(1),'facecolor',[.6 .6 .6],'edgecolor','k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-1 1.1],'xtick',[-1:0.5:1],'ylim',[0 25],'ytick',[0:12.5:25],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Performance Index')

subplot(3,2,3)
hold on
plot(data.summary.maxcorr(data.nsigind),data.crosspred.maxabove(data.nsigind),'o','markersize',6,...
  'markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6])
plot(data.summary.maxcorr(data.sigind),data.crosspred.maxabove(data.sigind),'x','markersize',6,...
  'markerfacecolor','k','markeredgecolor','k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[0 1],'xtick',[0:.5:1],'ylim',[0 1],'ytick',[0:.5:1],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Value of r in preferred direction')
ylabel('Value of r at +/- 45 degrees')

subplot(3,2,4)
hold on
plot(data.summary.maxcorr(data.nsigind),data.crosspred.maxbelow(data.nsigind),'o','markersize',6,...
  'markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6])
plot(data.summary.maxcorr(data.sigind),data.crosspred.maxbelow(data.sigind),'x','markersize',6,...
  'markerfacecolor','k','markeredgecolor','k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
'xlim',[0 1],'xtick',[0:.5:1],'ylim',[0 1],'ytick',[0:.5:1],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Value of r in preferred direction')

% figure describing how well the PD did, relative to how well it was
% possible to do

%indabv = (data.crosspred.maxabove - max(data.timecorr.abv')) ./ (data.crosspred.maxabove + max(data.timecorr.abv'));
%indblw = (data.crosspred.maxbelow - max(data.timecorr.blw')) ./ (data.crosspred.maxbelow + max(data.timecorr.blw'));
indabv = data.crosspred.maxabove./data.crosspred.actabove;
indblw = data.crosspred.maxbelow./data.crosspred.actbelow;

% a = hist(indabv(data.nsigind),[-1:0.2:1]);
% b = hist(indabv(data.unkind),[-1:0.2:1]);
% c = hist(indabv(data.classsigind),[-1:0.2:1]);

subplot(3,2,5)
hold on
hist(indabv(data.sigind),[0:.15:1.35]);
hist(indabv(data.nsigind),[0:.15:1.35]);
h = findobj(gca,'type','patch');
set(h(2),'facecolor','k','edgecolor',[.6 .6 .6])
set(h(1),'facecolor',[.6 .6 .6],'edgecolor','k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[0 1.5],'xtick',[0:0.5:1.5],'ylim',[0 50],'ytick',[0:25:50],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Performance above')
ylabel('Number of cells')

% a = hist(indblw(data.nsigind),[-1:0.2:1]);
% b = hist(indblw(data.unkind),[-1:0.2:1]);
% c = hist(indblw(data.classsigind),[-1:0.2:1]);
subplot(3,2,6)
hold on
hist(indblw(data.sigind),[0:.15:1.35]);
hist(indblw(data.nsigind),[0:.15:1.35]);
h = findobj(gca,'type','patch');
set(h(2),'facecolor','k','edgecolor',[.6 .6 .6])
set(h(1),'facecolor',[.6 .6 .6],'edgecolor','k')
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[0 1.5],'xtick',[0:0.5:1.5],'ylim',[0 50],'ytick',[0:25:50],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Performance below')
ylabel('Number of cells')

fie(fS4,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS4.xml')
print(fS4,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS4.eps')

%%%%%%%%%%%%%
% Figure S5 %
%%%%%%%%%%%%%

fS5 = figure;

subplot(3,2,1)
prefdirvec = 45:45:360;
circprefdirvec = 0:45:315;
eligible = data.pd(data.classsigind,:);
for i = 1:8
  dex = find(eligible(:,1) < prefdirvec(i) & eligible(:,1) >= circprefdirvec(i));
  meansx(i) = mean(eligible(dex,1));
  meansy(i) = mean(data.corr.rev(dex)); 
end

polar(deg2rad(data.pd(data.sigind,1)),data.summary.maxcorr(data.sigind)','kx')
hold on
a = polar(deg2rad(data.pd(data.nsigind,1)),data.summary.maxcorr(data.nsigind)','.');
set(a,'markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6);
a = polar(deg2rad([circprefdirvec circprefdirvec(1)]),[meansy meansy(1)]);
set(a,'linewidth',2,'color','k')
xlabel('preferred direction')
ylabel('r')

subplot(3,2,2)
hold on
plot(data.pd(data.sigind,2),data.summary.maxcorr(data.sigind),'kx','markersize',6)
plot(data.pd(data.nsigind,2),data.summary.maxcorr(data.nsigind),'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('tuning bandwidth')
ylabel('r')

subplot(3,2,3)
hold on
plot(data.fr(data.sigind),data.summary.maxcorr(data.sigind),'kx','markersize',6)
plot(data.fr(data.nsigind),data.summary.maxcorr(data.nsigind),'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('avg. firing rate')
ylabel('r')

subplot(3,2,4)
hold on
plot(data.numtrials(data.sigind),data.summary.maxcorr(data.sigind),'kx','markersize',6)
plot(data.numtrials(data.nsigind),data.summary.maxcorr(data.nsigind),'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('number of trials')
ylabel('r')

subplot(3,2,5)
hold on
plot(data.var.eye(data.sigind),data.summary.maxcorr(data.sigind),'kx','markersize',6)
plot(data.var.eye(data.nsigind),data.summary.maxcorr(data.nsigind),'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Variance in the eye movement')
ylabel('r')

subplot(3,2,6)
hold on
plot(data.var.sp(data.sigind),data.summary.maxcorr(data.sigind),'kx','markersize',6)
plot(data.var.sp(data.nsigind),data.summary.maxcorr(data.nsigind),'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6],'markersize',6)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Variance in the spike train')
ylabel('r')

fie(fS5,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS5.xml')
print(fS5,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureImport/FigureS5.eps')

%%%%%%%%%%%%%
% Figure S6 %
%%%%%%%%%%%%%

% A figure describing the PCA decomposition of the filters
filterbank = fliplr(data.filters.rev);

% flip the sign of the ones for 270 or 180
filterbank(find(data.prefdir == 270 | data.prefdir == 180),:) = -filterbank(find(data.prefdir == 270 | data.prefdir == 180),:);

% normalize the height
for i = 1:size(filterbank,1)
  filterbank(i,:) = filterbank(i,:)./norm(filterbank(i,:));
end

% pad the filters for the decomposition
filterbank = [zeros(141,200) filterbank zeros(141,200)];

% extract the appropriate filters
latency = data.filters.lat(find(~isnan(data.filters.lat)));
filterbank = filterbank(find(~isnan(data.filters.lat)),:);

% shift the filters
for i = 1:length(latency)
  filterbank(i,:) = circshift(filterbank(i,:),[1,-latency(i)]);
end

type = data.filters.type(find(~isnan(data.filters.lat)));
[u s v] = svd(filterbank);

PC1 = -filterbank*v(:,1);
PC2 = -filterbank*v(:,2);
PC3 = -filterbank*v(:,3);

fS6 = figure;
subplot(1,3,1)
hold on
plot(4:66,diag(s(4:end,4:end))./trace(s),'k.');
plot(1,s(1,1)./trace(s),'.','color',pretty('a'),'markersize',6)
plot(2,s(2,2)./trace(s),'.','color',pretty('i'),'markersize',6)
plot(3,s(3,3)./trace(s),'.','color',pretty('u'),'markersize',6)
hline(0.05)
set(gca,'xlim',[0 15],'ylim',[0 0.3])
ylabel('Variance described')

subplot(1,3,2)
hold on
hline(0)
plot(-550:550,-v(:,1),'color',pretty('a'),'linewidth',2)
plot(-550:550,-v(:,2),'color',pretty('i'),'linewidth',2)
plot(-550:550,v(:,3),'color',pretty('u'),'linewidth',2)
set(gca,'xlim',[-350 350],'ytick',[],'xtick',[-350 350])
xlabel('Time (ms)')

subplot(1,3,3)
hold on
plot(PC1(find(type == 1)),PC2(find(type == 1)),'.','color',pretty('i'))
plot(PC1(find(type == 2)),PC2(find(type == 2)),'.','color',pretty('a'))
plot(PC1(find(type == 3)),PC2(find(type == 3)),'.','color',pretty('u'))
set(gca,'box','off')
xlabel('PC1')
ylabel('PC2')

fie(fS6,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureS6.xml')
print(fS6,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureS6.eps')

%%%%%%%%%%%%%
% Figure S7 %
%%%%%%%%%%%%%
mu = [0.1 0.15 0.2];
sigma = [.01 .1];

% Uniform variance %

numtrials = 1000;
T = 1;
numunits = 250;  
  
for corrval = 1:3
 for run = 1:25;
  rho = mu(corrval)*ones(T,numtrials,numunits);
  nu = repmat(randn(1,numtrials),[T 1 numunits]);
  v = randn(T,numtrials,numunits);
  scale = ones(T,numtrials,numunits);
  
  M = squeeze(scale.*(sqrt(rho).*nu) + (sqrt(1-rho).*v));

  C = cov(M);
  NB = sum(C,2)./sqrt(diag(C)*sum(sum(C))); tmp = tril(C,-1); NN = tmp(find(tmp));
  [y x] = hist(NB,[-1:.01:1]); 
  static.NB(corrval,run,:) = y./sum(y);
  static.NN(corrval,run,:) = var(diag(C));
  [y x] = hist(diag(C),[.9:.05:1.1]);
  static.XX(corrval,run,:) = y./sum(y);
 end
end

colors = [pretty('a');pretty('i');pretty('u')];

fS7 = figure;
for i = 1:3,
  subplot(2,3,1)
  hold on
  errorbar(-1:.01:1,squeeze(mean(static.NB(i,:,:),2)),...
    squeeze(std(static.NB(i,:,:)./sqrt(25),[],2)),'color',colors(i,:));
  set(gca,'xlim',[-1 1],'ytick',[])
  
  subplot(2,3,4)
  hold on
  errorbar(.9:.05:1.1,squeeze(mean(static.XX(i,:,:),2)),...
    squeeze(std(static.XX(i,:,:)./sqrt(25),[],2)),'color',colors(i,:));
  set(gca,'xlim',[.85 1.15],'ytick',[])
  
end

% 25 neurons %

% k indexes sigma
% j indexes mu

for k = 1:2
  for j = 1:3
    for i = 1:1000,
      C = smartCov(25,mu(j),sigma(k),.01*j);
      NB = sum(C,2)./sqrt(diag(C)*sum(sum(C))); tmp = tril(C,-1); NN = tmp(find(tmp));
      [y x] = hist(NB,[-1:.01:1]); 
      r.NB(i,j,k,:) = y./sum(y); 
      r.NN(i,j,k) = var(NN); 
      [y x] = hist(diag(C),[0:.1:3]);
      r.m(i,j,k,:) = y./sum(y);
    end
  end
end

for k = 1:2
  for j = 1:3
    subplot(2,3,k+1)
    hold on
    errorbar([-1:.01:1],mean(squeeze(r.NB(:,j,k,:))),...
      std(squeeze(r.NB(:,j,k,:)))./sqrt(1000),'color',colors(j,:));
    set(gca,'xlim',[-1 1],'ytick',[]) 
    
   subplot(2,3,k+4)
   hold on
   errorbar([0:.1:3],mean(squeeze(r.m(:,j,k,:))),...
     std(squeeze(r.m(:,j,k,:)))./sqrt(1000),'color',colors(j,:));
    set(gca,'xlim',[0 3],'ytick',[])
  end
end

subplot(2,3,5), set(gca,'xlim',[0 1])

fie(fS7,'/Users/schoppik/Documents/Research/Papers/tbtv/FigureS7.xml')
print(fS7,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/FigureS7.eps')

%%%%%%%%%%%%%
% Figure S8 %
%%%%%%%%%%%%%

% The data to generate this figure were not saved.  Changing the
% analysis code will re-generate it.