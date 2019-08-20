fontsize = 14;
linewidth = 0.96;
ticklength = [.048 .048];

load example
load data

%%%%%%%%%%%%
% Figure 1 %
%%%%%%%%%%%%

windowsize = 10;

f1 = figure;
subplot(2,2,1)
hold on
x = -200:500;
y = [zeros(1,201) 20.*ones(1,500)];
plot(x,y,'color',pretty('k'),'linewidth',3,'linestyle','--')
eyex = [example.evel,nan(size(example.evel,1),1)]';
t = [repmat([-200:500 nan],1,size(example.evel,1))];
plot(t,eyex(:),'color',[.6 .6 .6])
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'tickdir','out',...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
plot(x,y,'color',pretty('k'),'linewidth',3,'linestyle','--') % in twice for legend
xlabel('Time (ms)')
ylabel('Eye Velocity (deg/s)')
h = legend('Target','Eye','location','northwest')
set(h,'FontSize',12);

subplot(2,2,2)
hold on
[x,y] = quickraster_phy(example.binaryspikes,200,500);
plot(x,y,'k')
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200,500],'xtick',[-200,0,250,500],'xticklabel',[-200 0 250 500],...
  'ylim',[1 size(example.binaryspikes,1)+1],'ytick',size(example.binaryspikes,1)+1,'yticklabel',size(example.binaryspikes,1),...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Time (ms)')
ylabel('Trial')

subplot(2,2,3)
hold on
eyex = [example.evel-repmat(mean(example.evel),size(example.evel,1),1),nan(size(example.evel,1),1)]';
t = [repmat([-200:500 nan],1,size(example.evel,1))];
plot(t,eyex(:),'color',[.6 .6 .6])
plot(-200:501,mean(eyex')+std(eyex'),'k','linewidth',2)
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'ylim',[-10 10],'ytick',[-10:10:10],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Time (ms)')
ylabel('Residual velocity (deg/s)')

subplot(2,2,4)
hold on
frx = [example.binaryspikes-repmat(mean(example.binaryspikes),size(example.binaryspikes,1),1),nan(size(example.binaryspikes,1),1)]';
t = [repmat([-200:500 nan],1,size(example.evel,1))];
plot(t,frx(:),'color',[.6 .6 .6])
plot(-200:501,mean(frx')+std(frx'),'k','linewidth',2)
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200 500],'xtick',[-200,0,250,500],'ylim',[-0.25 1],'ytick',[0:.5:1],...
  'fontsize',fontsize,'linewidth',linewidth,'TickLength',ticklength)
xlabel('Time (ms)')
ylabel('Residual spikes (imp/ms)')

fie(f1,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure1.xml')
print(f1,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/Figure1.eps')

%%%%%%%%%%%%
% Figure 2 %
%%%%%%%%%%%%

cutoff = .9;
f = [0:2047]*1000/2048;

for i = 1:141,
  maxcorr(i) = max(data.timecorr.rev(i,:));
  pottime = diff(find(data.timecorr.rev(i,:) >= (maxcorr(i)*cutoff)));
  maxcorrtime(i) = length(find(pottime == 1));
  maxcoh(i) = max(data.coherence(i,1:100));
  maxfreq(i).freq = f(find(data.coherence(i,1:100) >= (maxcoh(i)*cutoff)));
end


f2 = figure;
subplot(3,2,1)
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

subplot(3,2,2)
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

subplot(3,2,3)
hold on
errorbar(-200:500,mean(rand_tc),1.5*std(rand_tc),'color',[.6 .6 .6])
errorbar(-200:500,mean(tc),std(tc),'color','k')

set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-200,500],'xtick',[-200,0,250,500],'xticklabel',[-200 0 250 500],...
  'ylim',[-.5 1],'ytick',[-.5:.5:1])
xlabel('Time (ms)')
ylabel('Correlation coefficient')

subplot(3,2,4)
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


subplot(3,2,5)
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

subplot(3,2,6)
hold on
filtsize = log10(sum(abs(data.filters.rev')));
plot(filtsize(data.sigind),maxcorr(data.sigind),'k.')
plot(filtsize(data.nsigind),maxcorr(data.nsigind),'.','color',[.6 .6 .6])

set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'ylim',[0 1],'ytick',[0:.5:1],...
  'xlim',[1 3],'xtick',[1:1:3])
xlabel('log Filter size (deg/s)')
ylabel('Peak correlation')

fie(f2,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure2.xml')
print(f2,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/Figure2.eps')

%%%%%%%%%%%%
% Figure 3 %
%%%%%%%%%%%%

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

%%%%%%%%%%%%%
% now the NB correlations

maxcorr = max(data.timecorr.rev');
[jnk t] = max(data.timecorr.rev(sigind,:),[],2);
[jnk csiglag] = sort(t);

[jnk t] = max(data.timecorr.rev(nsigind,:),[],2);
[jnk cnsiglag] = sort(t);


f3 = figure;
subplot(2,3,1)
imagesc([psthnorm(sigind(siglag),:) ; nan(3,701);...
  psthnorm(nsigind(nsiglag),:)],[0 1])
g = colorbar;
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[1 701],'xtick',[1 201 450 701],'xticklabel',[-200 0 250 500],'ytick',[])
xlabel('Time (ms)')

subplot(2,3,2)
imagesc([data.timecorr.rev(sigind(csiglag),:) ; nan(3,701);...
  data.timecorr.rev(nsigind(cnsiglag),:)],[.25 .75])
g = colorbar;
set(gca,'tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[1 701],'xtick',[1 201 450 701],'xticklabel',[-200 0 250 500],'ytick',[])
xlabel('Time (ms)')

subplot(2,3,3)
hold on
plot(psth_vs_time(data.sigind),maxcorr(data.sigind),'k.')
plot(psth_vs_time(data.nsigind),maxcorr(data.nsigind),'.','color',[.6 .6 .6])
xlabel('PSTH vs. timecourse')
ylabel('Peak correlation')
set(gca,'PlotBoxAspectRatio',[1 1 1])

% Run list:
%
% 1: additive to both neurons, plot the firing rates
% 2: multiplicative on the first neuron, additive second
% 3: multiplicative first, multiplicative second
% 4: multiplicative first, additive second, MULTIPLICATIVE eye

% variables
ntrials = 500;
tlen = 701;
mnoisescale = 0.2;
anoisescale = 10;
numcells = 2;
output = zeros(ntrials,tlen,numcells);
Fhat = zeros(numcells,tlen);
color = [pretty('i'); pretty('u'); pretty('a')];

for run = 1:4
  for i = 1:numcells

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define filters and PSTHs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch i
      case 1
        rmean = sigmoid([100 0.05 400],[1:701]);
        %rmean = gaussian([350 20 100],[1:701]);
        F = gaussian([0 10 1],[-350:350]'); F = 0.5*F/norm(F);
      case 2
        rmean = gaussian([480 20 100],[1:701]);
        F = gaussian([0 10 1],[-350:350]'); F = 0.05*F/norm(F);
      case 3 % don't use this cell -- only for testing purposes
        rmean = sigmoid([100 0.05 400],[1:701]);
        F = gaussian([0 10 1],[-350:350]'); F = 0.5*F/norm(F);
    end

    if run == 1
      subplot(2,3,4); 
      hold on; 
      plot(-200:500,rmean,'-','color',color(i,:),'linewidth',2);       
      set(gca,'xlim',[-200 500],'ylim',[0 200],'xtick',[-200 0 250 500],'ytick',[0:50:200])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate the simulated neural activity %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % the mean firing rate
    base = repmat(rmean,ntrials,1);

    % multiplicative neural noise (noise that scales with the mean firing rate)
    mnn = mnoisescale.*base.*randn(ntrials,tlen); %corrnoise;

    % additive noise, added on to the neural noise
    ann = anoisescale.*randn(ntrials,tlen);

    % final neural activity
    switch run
      case 1 % additive noise to both neurons
        switch i
          case 1
            r = base + ann;
          case 2
            r = base + ann;
          case 3
            r = base + ann;
        end
      case {2} % multiplicative noise to the sigmoidal neuron
        switch i
          case 1
            r = base + mnn;
          case 2
            r = base + ann;
          case 3
            r = base + ann;
        end
        
      case 3 % multiplicative noise to both
        switch i
          case 1
            r = base + mnn;
          case 2
            r = base + 25.*mnn;
          case 3
            r = base + mnn;
        end
      case {4} % multiplicative noise to both
        switch i
          case 1
            r = base + mnn;
          case 2
            r = base + 50.*mnn;
          case 3
            r = base + ann;
        end
        

    end
    
    
    input(:,:,i) = r;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % derive the eye movements %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    command = myconv(r',repmat(F',ntrials,1)')';
 
    switch run
      case {1,2,3}
         % additive noise, added to the command
         acn = anoisescale.*randn(ntrials,tlen);
         e = command + acn;
      case 4
        % multiplicative noise (noise that scales with the command)
        mcn = mnoisescale.*command.*randn(ntrials,tlen);
        e = command + 2.5.*mcn;
    end

    output(:,:,i) = e;

  end


  for i = 1:numcells

    %%%%%%%%%%%%%%%%%%%%%%%
    % retrieve the filter %
    %%%%%%%%%%%%%%%%%%%%%%%

    % get the residuals
    varr(:,:,i) = input(:,:,i) - repmat(mean(input(:,:,i)),[ntrials 1]);
    
    S = mean(output,3);
    
    if run == 1
      subplot(2,3,4)
      hold on
      errorbar(-200:500,mean(S),std(S),'k')
    end
    
    vare = S - repmat(mean(S),[ntrials 1]);

    Fhat(i,:) = mean(myconv(varr(:,:,i)',vare'),2);
  
    %%%%%%%%%%%%%%%%%%%%%%%
    % make the prediction %
    %%%%%%%%%%%%%%%%%%%%%%%

    ehat = myconv(varr(:,:,i)',repmat(Fhat(i,:),ntrials,1)')';
  
    corrmat(i,:) = mycorr(ehat,S);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%
  % plot the predictions %
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  switch run
    case 1
      subplot(2,3,5)
      hold on
      plot(-200:500,corrmat(1,:),'color',[.6 .6 .6])
      set(gca,'xlim',[-200 500],'ylim',[-0.25 1],'xtick',[-200 0 250 500],'ytick',[0 .5 1])
    case 2
      subplot(2,3,5)
      hold on
      plot(-200:500,corrmat(1,:),'color',pretty('k'))
      set(gca,'xlim',[-200 500],'ylim',[-0.25 1],'xtick',[-200 0 250 500],'ytick',[0 .5 1])
    case 3
      subplot(2,3,6)
      hold on
      plot(-200:500,corrmat(1,:),'color',[.6 .6 .6])
      set(gca,'xlim',[-200 500],'ylim',[-0.25 1],'xtick',[-200 0 250 500],'ytick',[0 .5 1])
    case 4
      subplot(2,3,6)
      hold on
      plot(-200:500,corrmat(1,:),'-','color','k')
      set(gca,'xlim',[-200 500],'ylim',[-0.25 1],'xtick',[-200 0 250 500],'ytick',[0 .5 1])
  end
  

end % model runs


fie(f3,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure3.xml')
print(f3,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/Figure3.eps')

%%%%%%%%%%%%
% Figure 4 %
%%%%%%%%%%%%

load examplepair
sp1 = examplepair.spikes(find(examplepair.type == 180),:,1);
[x1 y1] = quickraster_phy(sp1(:,300:end),200,500);

sp2 = examplepair.spikes(find(examplepair.type == 180),:,2);
[x2 y2] = quickraster_phy(sp2(:,300:end),200,500);

prefdirvec = 45:45:360;
rundex = [1 2 3 4 5 8];

thresh = 1.5;

gooddir = 4;

XSAtime = squeeze(examplepair.rNB.AStime(rundex(gooddir),:,:));
XSBtime = squeeze(examplepair.rNB.BStime(rundex(gooddir),:,:));
randXSAtime = squeeze(examplepair.rNB.randAStime(rundex(gooddir),:,:));
randXSBtime = squeeze(examplepair.rNB.randBStime(rundex(gooddir),:,:));

XXAtime = squeeze(examplepair.rNB.AAtime(rundex(gooddir),:,:));
XXBtime = squeeze(examplepair.rNB.BBtime(rundex(gooddir),:,:));
NNtime = squeeze(examplepair.rNN.NNtime(rundex(gooddir),:,:));
NNtimecorr = squeeze(examplepair.rNN.NNtimecorr(rundex(gooddir),:,:));

tmp = mean(XSAtime) > thresh.*std(randXSAtime);
Adex{gooddir} = find(tmp);
  
tmp = mean(XSBtime) > thresh.*std(randXSBtime);
Bdex{gooddir} = find(tmp);

dex{gooddir} = intersect(Adex{gooddir},Bdex{gooddir});

if ~isempty(dex{gooddir})
  
  rXSAthresh = [mean(XSAtime(:,dex{gooddir}),2)'];
  rXSBthresh = [mean(XSBtime(:,dex{gooddir}),2)'];
  
  rXXAthresh = [mean(XXAtime(:,dex{gooddir}),2)'];
  rXXBthresh = [mean(XXBtime(:,dex{gooddir}),2)'];
  
  rNNthresh = [mean(NNtime(:,dex{gooddir}),2)'];
  
  rNNthreshcorr = mean(NNtimecorr(:,dex{gooddir}),2)';
  
  [jnk ted] = sort(mean(NNtimecorr),'descend');
  rNNthreshmaxcorr = mean(NNtimecorr(:,ted(1:length(dex{gooddir}))),2)';
end      

% get the average value of the correlation during the relevant time period
corr_during_overlap = mean(NNtimecorr(:,dex{4}).^2,2);
sorted = sort(NNtimecorr'.^2,'descend');
maximum_corr = mean(sorted(1:137,:));

f4 = figure;
% raster
subplot(2,3,1)
hold on
plot(x1,y1,'color',pretty('i'));
plot(x2,-y2,'color',pretty('u'));
set(gca,'xlim',[-200 500],'xtick',[-200 0 250 500],'ytick',[])
xlabel('Time (ms)')

% filters
subplot(2,3,2)
hold on
errorbar(-350:350,flipud(squeeze(mean(examplepair.rNB.randfB(4,:,:),2))),flipud(squeeze(std(examplepair.rNB.randfB(4,:,:),[],2))),'color',[.6 .6 .6]);
errorbar(-350:350,flipud(squeeze(mean(examplepair.rNB.fB(4,:,:),2))),flipud(squeeze(std(examplepair.rNB.fB(4,:,:),[],2))),'color',pretty('u'));
errorbar(-350:350,flipud(squeeze(mean(examplepair.rNB.fA(4,:,:),2))),flipud(squeeze(std(examplepair.rNB.fA(4,:,:),[],2))),'color',pretty('i'));
set(gca,'xlim',[-350 350],'xtick',[-350 0 350],'ylim',[-0.35 1.1],'ytick',[0:.5:1])
xlabel('Time from spike (ms)')
ylabel('Deg/s/spike')

% NB correlations
subplot(2,3,3)
hold on
errorbar(-200:500,squeeze(mean(examplepair.rNB.randBStime(4,:,:),2)),1.5.*squeeze(std(examplepair.rNB.randBStime(4,:,:),[],2)),'color',[.6 .6 .6])
errorbar(-200:500,squeeze(mean(examplepair.rNB.AStime(4,:,:),2)),squeeze(std(examplepair.rNB.AStime(4,:,:),[],2)),'color',pretty('i'))
errorbar(-200:500,squeeze(mean(examplepair.rNB.BStime(4,:,:),2)),squeeze(std(examplepair.rNB.BStime(4,:,:),[],2)),'color',pretty('u'))
rectangle('Position',[181 -0.35 137 1.35],'linewidth',2) 
plot(Adex{4}-200,.9*ones(1,length(Adex{4})),'color',pretty('i'),'linewidth',3)
plot(Bdex{4}-200,.8*ones(1,length(Bdex{4})),'color',pretty('u'),'linewidth',3)
set(gca,'xlim',[-200 500],'xtick',[-200 0 250 500],'ylim',[-0.35 1],'ytick',[0:.5:1])
ylabel('Correlation (r)')
xlabel('Time (ms)')

subplot(2,3,4)
hold on
xvals = [0.2:0.1:2.4];
[Ay,Ax] = hist(log10(XXAtime(:)),xvals);
[By,Bx] = hist(log10(XXBtime(:)),xvals);
[AyC,AxC] = hist(log10(rXXAthresh(:)),xvals);
[ByC,BxC] = hist(log10(rXXBthresh(:)),xvals);

bar(AxC,AyC./sum(AyC),'barwidth',.9,'facecolor','k','edgecolor',[.6 .6 .6])
bar(BxC,-ByC./sum(ByC),'barwidth',.9,'facecolor','k','edgecolor',[.6 .6 .6])

bar(Ax,Ay./sum(Ay),'barwidth',.9,'facecolor',pretty('i'))
bar(Bx,-By./sum(By),'barwidth',.9,'facecolor',pretty('u'))

set(gca,'xlim',[0 3],'xtick',[0:3],'ylim',[-.25 .25],'ytick',[])
xlabel('log Variation (deg/s)^2')
ylabel('Probability')

subplot(2,3,5)
hold on
xvals = [-20:2:30];
[y,x] = hist((NNtime(:)),xvals);
[yC,xC] = hist((rNNthresh(:)),xvals);

bar(xC,yC./sum(yC),'barwidth',.9,'facecolor','k','edgecolor',[.6 .6 .6])
bar(x,y./sum(y),'barwidth',.9,'edgecolor','k','facecolor',[.6 .6 .6])
set(gca,'xlim',[-25 35],'xtick',[-20:10:30],'ylim',[0 0.25],'ytick',[])
xlabel('Co-variation (deg/s)^2')

subplot(2,3,6)
hold on
plot([0 0.15],[0 0.15],'k--')
plot(corr_during_overlap,maximum_corr,'k.')
set(gca,'xlim',[0 0.15],'xtick',[0 0.15],'ylim',[0 0.15],'ytick',[0 0.15])
ylabel('Peak (r^2)')
xlabel('Overlap (r^2)')

fie(f4,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure4.xml')
print(f4,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/Figure4.eps')


%%%%%%%%%%%%
% Figure 5 %
%%%%%%%%%%%%

XX = pairdata.XXthresh;
XS = pairdata.XSthresh;
NN = pairdata.NNthresh;

% do some bootstrapping to figure out errors

EXX = fitlognorm(XX,1);
EXXhi = fitlognorm(XX,1) + sqrt(fitlognorm(XX,2));
%EXXlo = fitlognorm(XX,1) - sqrt(fitlognorm(XX,2));

ENN = median(NN);
NNs = sort(NN);
ENNhi = median(NN) + (NNs(round(length(NN)*3/4)) - NNs(round(length(NN)*1/4)));
%ENNlo = mean(NN) - std(NN);

EXShat = sqrt( (EXX + ENN.*[1:1000]) ./ (EXX .* [2:1001]));
EXShatHi = sqrt( (EXXhi + ENNhi.*[1:1000]) ./ (EXXhi .* [2:1001]));
%EXShatLo = sqrt( (EXXlo + ENNlo.*[1:1000]) ./ (EXXlo .* [2:1001]));


f5 = figure;

subplot(1,3,1)
[y,x] = hist(XS,100);

hold on
bar(x,(y./sum(y)./(x(2)-x(1))),'k');
plot([mean(XS) mean(XS)],[0 3.75],'color',[.6 .6 .6],'linewidth',2)
plot([mean(XS)+std(XS) mean(XS)+std(XS)],[0 3.75],':','color',[.6 .6 .6],'linewidth',2)
plot([mean(XS)-std(XS) mean(XS)-std(XS)],[0 3.75],':','color',[.6 .6 .6],'linewidth',2)
set(gca,'xlim',[-0.1 1],'ytick',[],'ylim',[0 3.75])
xlabel('Neuron-behavior correlation (r)')
ylabel('P')

subplot(1,3,2)
[y,x] = hist(log10(XX),100);

hold on
bar(x,(y./sum(y)./(x(2)-x(1))),'k');
plot(log10([fitlognorm(XX,1) fitlognorm(XX,1)]),[0 1.5],'color',[.6 .6 .6],'linewidth',2)
plot(log10([fitlognorm(XX,1) + sqrt(fitlognorm(XX,2)) fitlognorm(XX,1) + sqrt(fitlognorm(XX,2))]),[0 1.5],':','color',[.6 .6 .6],'linewidth',2)
set(gca,'xlim',[0 3],'ytick',[],'ylim',[0 1.5])
xlabel('log Variation (deg/s)^2')
ylabel('P')


subplot(1,3,3)
IQI = (NNs(round(length(NN)*3/4)) - NNs(round(length(NN)*1/4)));
[y,x] = hist(NN,300);

hold on
bar(x,(y./sum(y)./(x(2)-x(1))),'k');
plot([median(NN) median(NN)],[0 0.2],'color',[.6 .6 .6],'linewidth',2)
plot([median(NN)+IQI median(NN)+IQI],[0 0.2],':','color',[.6 .6 .6],'linewidth',2)
set(gca,'xlim',[-10 20],'ylim',[0 0.2],'ytick',[])
xlabel('Covariation (deg/s)^2')
ylabel('P')

fie(f5,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure5.xml')
print(f5,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/Figure5.eps')

%%%%%%%%%%%%
% Figure 6 %
%%%%%%%%%%%%

f6 = figure ;
hold on
errorbar(2:1001,mean(XS)*ones(1,1000),std(XS)*ones(1,1000),'k')
errorbar(2:1001,EXShat,EXShatHi-EXShat,'color',[.6 .6 .6])
set(gca,'xlim',[2 10e2],'xtick',[10e0 10e1 10e2],'xscale','log','ylim',[0 1])
xlabel('Number of neurons')
ylabel('Est. neuron-behavior correlation (r)')

fie(f6,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure6.xml')
print(f6,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/Figure6.eps')


%%%%%%%%%%%%
% Figure 7 %
%%%%%%%%%%%%

% A figure describing the filters
filterbank = fliplr(data.filters.rev);
% flip the sign of the ones for 270 or 180
filterbank(find(data.prefdir == 270 | data.prefdir == 180),:) = -filterbank(find(data.prefdir == 270 | data.prefdir == 180),:);

% normalize the height
for i = 1:size(filterbank,1)
  filterbank(i,:) = filterbank(i,:)./norm(filterbank(i,:));
end

f7 = figure;
subplot(2,3,1)
hold on
plot([-350 350],[0 0],'k:')
plot([0 0],[-0.15 .15],'k:')
ind = find(data.filters.type == 1);
for i = 1:length(ind)
  x = [-350:350] - data.filters.lat(ind(i));
  y = filterbank(ind(i),:);
  plot(x,y,'color',[.6 .6 .6])
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
  plot(x,y,'color',[.6 .6 .6])
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
  plot(x,y,'color',[.6 .6 .6])
end
set(gca,'box','off','tickdir','out','PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','auto',...
  'xlim',[-350 350],'xtick',[-350 0 350],'ylim',[-0.15 .15],'ytick',[])

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

subplot(2,3,5)
hold on
plot([-.1 1],[-.1 1],'k--')
plot(pairdata.rSC.prs.m(sigdex),pairdata.NNcorr.m(sigdex),'k.')
set(gca,'xlim',[-.1 1],'ylim',[-.1 1])
xlabel('Pursuit (r)')
ylabel('Filtered spikes (r)')

xvals = [0:3];
ed = pairdata.electrodedist(sigdex);
for i = 1:4
  
  m(i) = mean(pairdata.NNcorr.m(sigdex(find(ed == i-1))));
  s(i) = std(pairdata.NNcorr.m(sigdex(find(ed == i-1))));
end

subplot(2,3,6)
hold on
plot(pairdata.electrodedist(sigdex),pairdata.NNcorr.m(sigdex),'.','color',[.6 .6 .6])
errorbar(xvals,m,s,'ko')
set(gca,'xlim',[-1 4],'xtick',[0:1:3],'ylim',[-0.05 0.45],'ytick',[0:.1:.4])
xlabel('Electrode spacing')
ylabel('Filtered spike train correlation')

fie(f7,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure7.xml')
print(f7,'-depsc','/Users/schoppik/Documents/Research/Papers/tbtv/Figure7.eps')

%%%%%%%%%%%%
% Figure 8 %
%%%%%%%%%%%%
% temporary, before reshuffling

corrMean = mean(pairdata.thresh.NNcorrsum);
prodMean = mean(pairdata.thresh.prodsum);

theoreticalWidth = 2*mean(tmps);
actualWidth = 2*std(pairdata.thresh.prodsum);

height = 2*std(pairdata.thresh.NNcorrsum);


f8 = figure;
subplot(1,2,1)
hold on

plot([-.1 .5],[-.1 .5],':','color',[.6 .6 .6])
plot(pairdata.thresh.prodsum,pairdata.thresh.NNcorrsum,'k.')

for i = 0:4
  tmpM(i + 1) = mean(pairdata.thresh.NNcorrsum(find(pairdata.electrodedist == i)) - ...
    pairdata.thresh.prodsum(find(pairdata.electrodedist == i)));
  tmpS(i + 1) = std(pairdata.thresh.NNcorrsum(find(pairdata.electrodedist == i)) - ...
    pairdata.thresh.prodsum(find(pairdata.electrodedist == i)));
end

% theoretical fit
rectangle('curvature',[1 1],'position',[corrMean - (theoreticalWidth/2) corrMean - (height/2) ...
  theoreticalWidth height],'edgecolor','r','linewidth',2);

% actual fit
rectangle('curvature',[1 1],'position',[prodMean - (actualWidth/2) corrMean - (height/2) ...
  theoreticalWidth height],'edgecolor','k','linewidth',2);


set(gca,'plotboxaspectratiomode','manual','box','off','xlim',[-.1 .5],'ylim',[-.1 .5],...
  'xtick',[0:.25:.5],'ytick',[0:.25:.5])
xlabel('rNB1 * rNB2');
ylabel('rNN')

subplot(1,2,2)
hold on
plot([-1 4],[0 0],':','color',[.6 .6 .6]);
plot(pairdata.electrodedist,pairdata.thresh.NNcorrsum-pairdata.thresh.prodsum,'.','color',[.6 .6 .6])
errorbar(0:4,tmpM,tmpS,'o','color','k')
set(gca,'plotboxaspectratiomode','manual','box','off','xlim',[-1 4],'ylim',[-.3 .3],...
  'ytick',[-.3:.1:.3],'xtick',[0:1:3])
xlabel('Inter-electrode spacing')
ylabel('rNN - (rNB1 * rNB2)')

fie(f8,'/Users/schoppik/Documents/Research/Papers/tbtv/Figure8.xml')


