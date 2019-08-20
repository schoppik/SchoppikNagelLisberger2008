%
%
%
%          by: david schoppik
%        date: 8/1/2007
%     purpose: to generate figures

load example
load data
load examplepair
load pairdata
load pair098

% figure 1 %

windowsize = 10;

f1 = figure;
subplot(2,2,1)
hold on
x = -200:500;
y = [zeros(1,201) 20.*ones(1,500)];
plot(x,y,'color',pretty('k'),'linewidth',3,'linestyle',':')
eyex = [example.evel,nan(size(example.evel,1),1)]';
t = [repmat([-200:500 nan],1,size(example.evel,1))];
plot(t,eyex(:),'color',[.6 .6 .6])
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'tickdir','out')
plot(x,y,'color',pretty('k'),'linewidth',3,'linestyle','--') % in twice for legend
xlabel('Time (ms)')
ylabel('Eye Velocity (deg/s)')

subplot(2,2,2)
hold on
[x,y] = quickraster_phy(example.binaryspikes,200,500);
plot(x,y,'k')
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200,500],'xtick',[-200,0,250,500],'xticklabel',[-200 0 250 500],...
  'ylim',[1 size(example.binaryspikes,1)+1],'ytick',size(example.binaryspikes,1)+1,'yticklabel',size(example.binaryspikes,1))
xlabel('Time (ms)')
ylabel('Trial')

subplot(2,2,3)
hold on

eyex = [example.evel-repmat(mean(example.evel),size(example.evel,1),1),nan(size(example.evel,1),1)]';
t = [repmat([-200:500 nan],1,size(example.evel,1))];
plot(t,eyex(:),'color',[.6 .6 .6])
plot(-200:501,mean(eyex')+std(eyex'),'k','linewidth',2)
plot(-200:501,mean(eyex')-std(eyex'),'k','linewidth',2)
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'ylim',[-10 10],'ytick',[-10:10:10])
xlabel('Time (ms)')
ylabel('Residual velocity (deg/s)')

subplot(2,2,4)
hold on
frx = [example.binaryspikes-repmat(mean(example.binaryspikes),size(example.binaryspikes,1),1),nan(size(example.binaryspikes,1),1)]';
t = [repmat([-200:500 nan],1,size(example.evel,1))];
plot(t,frx(:),'color',[.6 .6 .6])
plot(-200:501,mean(frx')+std(frx'),'k','linewidth',2)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'ylim',[-0.25 1],'ytick',[0:.5:1])
xlabel('Time (ms)')
ylabel('Residual spikes (imp/ms)')

% PSTH heatmap
PSTH = windowspikedata(example.binaryspikes,10);
f1a = figure;
imagesc(sum(PSTH)./max(sum(PSTH)));

% figure 2 %

f2 = figure;
subplot(2,2,1)
imagesc(example.fefr,[-.7 .7])
hold on
plot([1 700],[1 700],'color','k','linewidth',2)
g = colorbar;
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'ylim',[1 700],'ytick',[1 200 450 700],'yticklabel',[-200 0 250 500],...
  'xlim',[1 700],'xtick',[1 200 450 700],'xticklabel',[-200 0 250 500],'xaxislocation','top')
xlabel('Firing rate')
ylabel('Eye velocity')
h= get(gcf,'children');
set(h(1),'box','off','tickdir','out');
set(g,'ytick',[-.7:.35:.7])

subplot(2,2,2)
hold on
errorbar(-350:350,flipud(mean(example.filterbank.randvar,2)),...
  flipud(std(example.filterbank.randvar,0,2)),...
  flipud(std(example.filterbank.randvar,0,2)),'color',[.6 .6 .6]);

errorbar(-350:350,flipud(mean(example.filterbank.var,2)),...
  flipud(std(example.filterbank.var,0,2)),...
  flipud(std(example.filterbank.var,0,2)),'color','k');

set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-350 350],'xtick',[-350:175:350],'ylim',[-.25 1],'ytick',[-.25:.25:1],...
  'xaxislocation','top','yaxislocation','right')
hline(0)
vline(0)

xlabel('Time (ms)')
ylabel('Deg/s per spike')

rand_tc = example.test.shuftimecorr;
tc = example.test.timecorr;

subplot(2,2,3)
hold on
errorbar(-200:500,mean(rand_tc),1.5*std(rand_tc),'color',[.6 .6 .6])
errorbar(-200:500,mean(tc),std(tc),'color','k')

set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200,500],'xtick',[-200,0,250,500],'xticklabel',[-200 0 250 500],...
  'ylim',[-.5 1],'ytick',[-.5:.5:1])
xlabel('Time (ms)')
ylabel('Correlation coefficient')

subplot(2,2,4)
nfft = 2048;
A = (example.evel - repmat(mean(example.evel),[size(example.evel,1) 1]))';
B = (example.binaryspikes - repmat(mean(example.binaryspikes),[size(example.binaryspikes,1) 1]))';
fAB = mean(fft(A,nfft).*conj(fft(B,nfft)),2);
fABs = mean(fft(A,nfft).*conj(fft(B(:,randperm(size(B,2))),nfft)),2);
fAA = mean(fft(A,nfft).*conj(fft(A,nfft)),2);
fBB = mean(fft(B,nfft).*conj(fft(B,nfft)),2);
coh = (abs(fAB).^2)./(fAA.*fBB);
scoh = (abs(fABs).^2)./(fAA.*fBB);
f = [0:nfft-1]*1000/nfft;  

hold on
bar(f(1:34),coh(1:34),'facecolor','k','edgecolor','k','linewidth',2)
bar(f(1:34),scoh(1:34),'facecolor',[.6 .6 .6],'edgecolor',[.6 .6 .6],'linewidth',2)
set(gca,'box','off','xtick',[0:2:16],'xlim',[0 16],'ylim',[0 0.3],'ytick',[0 0.3])
xlabel('Frequency (Hz)')
ylabel('Coherence')

axes('position',get(gca,'position'))
plot(f(1:34),fAA(1:34)./sum(fAA),'r','linewidth',3)
set(gca,'color','none','xtick',[],'xlim',[0 16],'ylim',[0 0.1],'ytick',[],'yaxislocation','right')
ylabel('Power in eye')

% figure 3 %

% smooth and normalize the PSTHs
f = (normalizedgaussian([0 25 1],-100:100));

for i = 1:141, 
  tmp = conv(f,(data.psth(i,:)./max(data.psth(i,:)))); 
  Fpsth(i,:) = tmp(101:end-100); 
end

psth_vs_time = mycorr(Fpsth',data.timecorr.rev');

% now, prepare to show the spikes
for i = 1:141
  psth(i,:) = mean(windowspikedata(data.binaryspikes{i},10));
end

psthnorm = (psth./repmat(max(psth,[],2),[1 701]));
sigind = data.sigind;
nsigind = data.nsigind;

% sort by the time to peak
[jnk t] = max(psthnorm(sigind,:),[],2);
[jnk siglag] = sort(t);

[jnk t] = max(psthnorm(nsigind,:),[],2);
[jnk nsiglag] = sort(t);

% now the NB correlations
maxcorr = max(data.timecorr.rev');
[jnk t] = max(data.timecorr.rev(sigind,:),[],2);
[jnk csiglag] = sort(t);

[jnk t] = max(data.timecorr.rev(nsigind,:),[],2);
[jnk cnsiglag] = sort(t);


f3 = figure;

subplot(1,3,1)
hold on
hist(maxcorr(data.sigind),[0:.075:825]);
hist(maxcorr(data.nsigind),[0:.075:.525]);
h = findobj(gca,'type','patch');
set(h(2),'facecolor','k','edgecolor',[.6 .6 .6])
set(h(1),'facecolor',[.6 .6 .6],'edgecolor','k')
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'ylim',[0 30],'ytick',[0:15:30],...
  'xlim',[0 .8],'xtick',[0:.4:.8])
plot([mean(maxcorr(data.sigind)) mean(maxcorr(data.sigind))],[0 30],'k:','linewidth',2);
plot([mean(maxcorr(data.nsigind)) mean(maxcorr(data.nsigind))],[0 30],':','color',[.6 .6 .6],'linewidth',2);
 xlabel('Peak correlation')
ylabel('Units')

subplot(1,3,2)
imagesc([psthnorm(sigind(siglag),:) ; nan(3,701);...
  psthnorm(nsigind(nsiglag),:)],[0 1])
g = colorbar;
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[1 701],'xtick',[1 201 450 701],'xticklabel',[-200 0 250 500],'ytick',[])
xlabel('Time (ms)')

subplot(1,3,3)
imagesc([data.timecorr.rev(sigind(csiglag),:) ; nan(3,701);...
  data.timecorr.rev(nsigind(cnsiglag),:)],[.25 .75])
g = colorbar;
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[1 701],'xtick',[1 201 450 701],'xticklabel',[-200 0 250 500],'ytick',[])
xlabel('Time (ms)')

% figure 4

% A figure describing the filters
filterbank = fliplr(data.filters.rev);
% flip the sign of the ones for 270 or 180
filterbank(find(data.prefdir == 270 | data.prefdir == 180),:) = -filterbank(find(data.prefdir == 270 | data.prefdir == 180),:);

% normalize the height
for i = 1:size(filterbank,1)
  filterbank(i,:) = filterbank(i,:)./norm(filterbank(i,:));
end

f4 = figure;
subplot(2,3,1)
hold on
plot([-350 350],[0 0],'k:')
plot([0 0],[-0.15 .15],'k:')
ind = find(data.filters.type == 1);
for i = 1:length(ind)
  x = [-350:350] - data.filters.lat(ind(i));
  y = filterbank(ind(i),:);
  plot(x,y,'color',pretty('i'),'linewidth',.25)
end

set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-350 350],'xtick',[-350 0 350],'ylim',[-0.15 .15],'ytick',[])


subplot(2,3,2)
hold on
plot([-350 350],[0 0],'k:')
plot([0 0],[-0.15 .15],'k:')
ind = find(data.filters.type == 2);
for i = 1:length(ind)
  x = [-350:350] - data.filters.lat(ind(i));
  y = filterbank(ind(i),:);
  plot(x,y,'color',pretty('a'),'linewidth',.25)
end

set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-350 350],'xtick',[-350 0 350],'ylim',[-0.15 .15],'ytick',[])
xlabel('Lag (ms)')

subplot(2,3,3)
hold on
plot([-350 350],[0 0],'k:')
plot([0 0],[-0.15 .15],'k:')
ind = find(data.filters.type == 3);
for i = 1:length(ind)
  x = [-350:350] - data.filters.lat(ind(i));
  y = filterbank(ind(i),:);
  plot(x,y,'color',pretty('u'),'linewidth',.25)
end
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-350 350],'xtick',[-350 0 350],'ylim',[-0.15 .15],'ytick',[])

% the "stacked" bar chart has to be created in Phyplot.  
subplot(2,3,4)
hold on
hist(data.filters.lat(find(data.filters.type == 1)),[-160:20:160]);
hist(data.filters.lat(find(data.filters.type == 2)),[-160:20:160]);
hist(data.filters.lat(find(data.filters.type == 3)),[-160:20:160]);
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 200],'xtick',[-200:50:200],'ylim',[0 15],'ytick',[0:5:15])
xlabel('Moment of Maximal Influence (ms)')
ylabel('Filters')

sigdex = find(pairdata.rSC.prs.sig);
bothsigdex = find(pairdata.rSC.prs.sig & pairdata.rSC.fix.sig);

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

subplot(2,3,5)
hold on
hline(0)
plot(-550:550,-v(:,1),'color',pretty('a'),'linewidth',2)
plot(-550:550,-v(:,2),'color',pretty('i'),'linewidth',2)
plot(-550:550,v(:,3),'color',pretty('u'),'linewidth',2)
set(gca,'xlim',[-350 350],'ytick',[],'xtick',[-350 350])
xlabel('Time (ms)')

subplot(2,3,6)
hold on
plot(PC1(find(type == 1)),PC2(find(type == 1)),'.','color',pretty('i'))
plot(PC1(find(type == 2)),PC2(find(type == 2)),'.','color',pretty('a'))
plot(PC1(find(type == 3)),PC2(find(type == 3)),'.','color',pretty('u'))
set(gca,'box','off')
xlabel('PC1')
ylabel('PC2')

% figure 5
sp1 = examplepair.spikes(find(examplepair.type == 180),:,1);
[x1 y1] = quickraster_phy(sp1(:,300:end),200,500);

sp2 = examplepair.spikes(find(examplepair.type == 180),:,2);
[x2 y2] = quickraster_phy(sp2(:,300:end),200,500);

prefdirvec = 45:45:360;
rundex = [1 2 3 4 5 8];
thresh = 1.5;

% extract the relevant values at the appropriate times
gooddir = 4;
NNtimecorr = squeeze(examplepair.rNN.NNtimecorr(rundex(gooddir),:,:));
NBA = squeeze(examplepair.rNB.AStime(rundex(gooddir),:,:));
NBB = squeeze(examplepair.rNB.BStime(rundex(gooddir),:,:));

XSAtime = squeeze(examplepair.rNB.AStime(rundex(gooddir),:,:));
XSBtime = squeeze(examplepair.rNB.BStime(rundex(gooddir),:,:));
randXSAtime = squeeze(examplepair.rNB.randAStime(rundex(gooddir),:,:));
randXSBtime = squeeze(examplepair.rNB.randBStime(rundex(gooddir),:,:));

tmp = mean(XSAtime) > thresh.*std(randXSAtime);
Adex{gooddir} = find(tmp);
tmp = mean(XSBtime) > thresh.*std(randXSBtime);
Bdex{gooddir} = find(tmp);
dex{gooddir} = intersect(Adex{gooddir},Bdex{gooddir});
rNN = mean(NNtimecorr(:,dex{gooddir}),1)';
rNBstar = mean(NBA(:,dex{gooddir}),1)' .* mean(NBB(:,dex{gooddir}),1)';

f5 = figure;
% raster
subplot(2,3,1)
hold on
plot(x1,y1,'color',pretty('i'));
plot(x2,-y2,'color',pretty('u'));
set(gca,'xlim',[-200 500],'xtick',[-200 0 250 500],'ytick',[])
xlabel('Time (ms)')

% NB correlations
subplot(2,3,2)
hold on
errorbar(-200:500,squeeze(mean(examplepair.rNB.AStime(4,:,:),2)),squeeze(std(examplepair.rNB.AStime(4,:,:),[],2)),'color',pretty('i'))
errorbar(-200:500,squeeze(mean(examplepair.rNB.BStime(4,:,:),2)),squeeze(std(examplepair.rNB.BStime(4,:,:),[],2)),'color',pretty('u'))
errorbar(-200:500,squeeze(mean(examplepair.rNN.NNtimecorr(4,:,:),2)),squeeze(std(examplepair.rNN.NNtimecorr(4,:,:),[],2)),'color','k')
plot(Adex{4}-200,.9*ones(1,length(Adex{4})),'color',pretty('i'),'linewidth',3)
plot(Bdex{4}-200,.8*ones(1,length(Bdex{4})),'color',pretty('u'),'linewidth',3)
set(gca,'xlim',[-200 500],'xtick',[-200 0 250 500],'ylim',[-0.35 1],'ytick',[0:.5:1])
ylabel('Correlation (r)')
xlabel('Time (ms)')

% rNB* versus rNN
subplot(2,3,3)
hold on
plot([-.1 .4],[-.1 .4],'k:')
plot(rNBstar,rNN,'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6])
plot(mean(rNBstar),mean(rNN),'k*')
set(gca,'xlim',[-.1 .4],'ylim',[-.1 .4])
xlabel('rNB*')
ylabel('rNN')
% now, the second pair

sp1 = pair098.spikes(find(pair098.type == 315),:,1);
[x1 y1] = quickraster_phy(sp1(:,300:end),200,500);

sp2 = pair098.spikes(find(pair098.type == 315),:,2);
[x2 y2] = quickraster_phy(sp2(:,300:end),200,500);

prefdirvec = 45:45:360;
rundex = [1 7 8];
thresh = 1.5;

% extract the relevant values at the appropriate times
gooddir = 2;
NNtimecorr = squeeze(pair098.rNN.NNtimecorr(rundex(gooddir),:,:));
NBA = squeeze(pair098.rNB.AStime(rundex(gooddir),:,:));
NBB = squeeze(pair098.rNB.BStime(rundex(gooddir),:,:));

XSAtime = squeeze(examplepair.rNB.AStime(rundex(gooddir),:,:));
XSBtime = squeeze(examplepair.rNB.BStime(rundex(gooddir),:,:));
randXSAtime = squeeze(examplepair.rNB.randAStime(rundex(gooddir),:,:));
randXSBtime = squeeze(examplepair.rNB.randBStime(rundex(gooddir),:,:));

tmp = mean(XSAtime) > thresh.*std(randXSAtime);
Adex{gooddir} = find(tmp);
tmp = mean(XSBtime) > thresh.*std(randXSBtime);
Bdex{gooddir} = find(tmp);
dex{gooddir} = intersect(Adex{gooddir},Bdex{gooddir});
rNN = mean(NNtimecorr(:,dex{gooddir}),1)';
rNBstar = mean(NBA(:,dex{gooddir}),1)' .* mean(NBB(:,dex{gooddir}),1)';

% raster
subplot(2,3,4)
hold on
plot(x1,y1,'color',pretty('i'));
plot(x2,-y2,'color',pretty('u'));
set(gca,'xlim',[-200 500],'xtick',[-200 0 250 500],'ytick',[])
xlabel('Time (ms)')

% NB correlations
subplot(2,3,5)
hold on
errorbar(-200:500,squeeze(mean(pair098.rNB.AStime(7,:,:),2)),squeeze(std(pair098.rNB.AStime(7,:,:),[],2)),'color',pretty('i'))
errorbar(-200:500,squeeze(mean(pair098.rNB.BStime(7,:,:),2)),squeeze(std(pair098.rNB.BStime(7,:,:),[],2)),'color',pretty('u'))
plot(-200:500,squeeze(mean(pair098.rNN.NNtimecorr(7,:,:),2)),'color','k','linewidth',2)
plot(Adex{2}-200,.9*ones(1,length(Adex{2})),'color',pretty('i'),'linewidth',3)
plot(Bdex{2}-200,.8*ones(1,length(Bdex{2})),'color',pretty('u'),'linewidth',3)
set(gca,'xlim',[-200 500],'xtick',[-200 0 250 500],'ylim',[-0.35 1],'ytick',[0:.5:1])
ylabel('Correlation (r)')
xlabel('Time (ms)')

% rNB* versus rNN
subplot(2,3,6)
hold on
plot([-.1 .4],[-.1 .4],'k:')
plot(rNBstar,rNN,'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6])
plot(mean(rNBstar),mean(rNN),'k*')
set(gca,'xlim',[-.1 .4],'ylim',[-.1 .4])
xlabel('rNB*')
ylabel('rNN')

% figure 6
[y x] = hist(pairdata.thresh.XX,15);
dist = x(2)-x(1);
params = lognfit(pairdata.thresh.XX);
fitx = 0:.5:185;
fity = lognpdf(fitx,params(1),params(2));
fity = fity./sum(fity).*dist;

f6 = figure;
subplot(1,3,1)
hold on
bar(x,y./sum(y),'barwidth',.9,'facecolor',[.6 .6 .6],'edgecolor','k')
plot(fitx,fity,'k','linewidth',3)
xlabel('sigmaN')
ylabel('Probability')

subplot(1,3,2)
hold on
plot([-.1 .5],[-.1 .5],'k:')
plot(pairdata.thresh.prodsum,pairdata.thresh.NNcorrsum,'.',...
  'markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6]);
% plot lines to guide the ellipse
mx = mean(pairdata.thresh.prodsum);
my = mean(pairdata.thresh.NNcorrsum);
sx = 2.*std(pairdata.thresh.prodsum);
sy = 2.*std(pairdata.thresh.NNcorrsum);
plot([mx mx],[my - sy my + sy],'k')
plot([mx-sx mx+sx],[my my],'k')

set(gca,'xlim',[-.1 .5],'ylim',[-.1 .5])
xlabel('rNB*')
ylabel('rNN')

subplot(1,3,3)
hold on
plot([-1 4],[0 0],'k:')
m = [];
s = [];
for i = 0:3
  dex = find(pairdata.electrodedist == i);
  plot(i,pairdata.thresh.NNcorrsum(dex) - pairdata.thresh.prodsum(dex),'.','markerfacecolor',[.6 .6 .6],'markeredgecolor',[.6 .6 .6]);
  m(i+1) = mean(pairdata.thresh.NNcorrsum(dex) - pairdata.thresh.prodsum(dex));
  s(i+1) = std(pairdata.thresh.NNcorrsum(dex) - pairdata.thresh.prodsum(dex));
end

errorbar(0:3,m,s,'k.');
set(gca,'xlim',[-1 4],'ylim',[-.3 .3])
xlabel('Electrode spacing')
ylabel('rNN - rNB*')

% figure 7
% run a simulation
% change how thoroughly we want to map the space
numtrials = 2500;
numRunsToAverage = 50;
numBins = 50;

% real rho = 0.1835, real var = 0.0125
desiredRho = 0.10705;
desiredRhoVar = 0.014;
threshold = 0.11;

% fit the data
params = lognfit(pairdata.thresh.XX);
lognormMu = params(1);
lognormSigma = params(2);

numunitsVec = round(logspace(1,2,numBins));
noiseLevelVec = [0 logspace(-2.5,-1,numBins)];

rNBMat = zeros(length(noiseLevelVec),length(numunitsVec));

noiseInd = 1;
for noiseLevel = noiseLevelVec
  unitInd = 1;
  for numunits = numunitsVec
   
    % clear out the old junk
    tmpm = zeros(1,numRunsToAverage);
    for run = 1:numRunsToAverage
          
      corrMat = smartCov(numunits,desiredRho,desiredRhoVar,threshold);
      randomDraws = mvnrnd(zeros(numtrials,numunits),corrMat);
      
      % now, turn the matrix of correlations into a covariance matrix
      randomWeights = sqrt(lognrnd(lognormMu,lognormSigma,[1 numunits]));
      
      % scale randomDraws by the variance
      randomDraws = repmat(randomWeights,[numtrials 1]).*randomDraws;
             
      actualSum = sum(randomDraws,2);
      
      % noise scale should be a fraction of the variance in the actual sum
      % a noise level of 1 should be the same size as the variance in the
      % actual sum
      noiseScale = noiseLevel * std(actualSum);
      noisyBehavior = actualSum + (actualSum .* (noiseScale.*randn(size(actualSum))));

      % find out the pairwise values
      allCorrelations = mycorr(repmat(sum(noisyBehavior,2),[1 numunits]),randomDraws);
      productofCorrelations = tril(allCorrelations*allCorrelations',-1);
      tmpm(run) = median(productofCorrelations(find(productofCorrelations)));
    end
    
    
    rNBMat(noiseInd,unitInd) = mean(tmpm);
    unitInd = unitInd + 1;

  end
  noiseInd = noiseInd + 1;
  disp(num2str((noiseInd-1)/length(noiseLevelVec)))
end


