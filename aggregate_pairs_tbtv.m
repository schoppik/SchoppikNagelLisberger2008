
% function to aggregate the pairs, post-analysis, into a data structure
% assumes that the pair files are in the path

% significance threshold
thresh = 1.5;
minnumtimepts = 20;

prefdirvec = 45:45:360;

% preallocation
pairdata = [];

pairdata.rSC.fix.m = [];
pairdata.rSC.prs.m = [];
pairdata.rSC.fix.s = [];
pairdata.rSC.prs.s = [];
pairdata.rSC.fix.sig = [];
pairdata.rSC.prs.sig = [];

pairdata.pd = [];
pairdata.ind = [];
pairdata.numdirs = [];
pairdata.electrodedist = [];
pairdata.analid = [];
pairdata.timepts = [];

pairdata.filtdp.m = [];
pairdata.filtdp.s = [];

pairdata.XSthresh = [];
pairdata.XXthresh = [];
pairdata.NNthresh = [];
pairdata.thresh.XXA = [];
pairdata.thresh.XXB = [];
pairdata.thresh.NNcorr = [];
pairdata.thresh.randNNcorr = [];
pairdata.thresh.maxcorr = [];
pairdata.thresh.frac = [];
pairdata.thresh.XX = [];


pairdata.thresh.prod = [];

pairdata.thresh.prodsum = [];
pairdata.thresh.NNcorrsum = [];

pairdata.uNNthresh = [];

pairdata.NNcorr.m = [];
pairdata.NNcorr.s = [];


% begin the pair loop
for pair_ind = [1:104];

  % load the unit and place it into a placeholder
  load(sprintf('pair%.3d',pair_ind));
  eval(sprintf('pair = pair%.3d;',pair_ind));

  pairdata.pd = [pairdata.pd abs(pair.prefdir(1)-pair.prefdir(2))];
  
  randex = find(pair.run);
  issig = [];
  
  if ~isempty(randex) % did we run any at all?
    
      % check for significance, by direction
    for i = 1:length(randex)
      FfA = sum(abs(fft(squeeze(pair.rNB.fA(randex(i),:,:)).^2)),2);
      FrandfA = sum(abs(fft(squeeze(pair.rNB.randfA(randex(i),:,:)).^2)),2);
      
      FfB = sum(abs(fft(squeeze(pair.rNB.fB(randex(i),:,:)).^2)),2);
      FrandfB = sum(abs(fft(squeeze(pair.rNB.randfB(randex(i),:,:)).^2)),2);
      
      issig(i) = and(ttest(FfA,FrandfA,.05,'right'),ttest(FfB,FrandfB,.05,'right'));
    end

    rundex = randex(find(issig));
    nrundex = randex(find(~issig)); 
    
    if ~isempty(rundex) % were there any significant filters

      %%%%%%%%%%%%%%%%%%%%%
      % housekeeping info %
      %%%%%%%%%%%%%%%%%%%%%
      pairdata.ind = [pairdata.ind pair_ind];
      pairdata.numdirs = [pairdata.numdirs length(rundex)];
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Threshold values in time %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      filtdpm = [];
      filtdps = [];
      NNm = [];
      NNs = [];
      
      for gooddir = 1:length(rundex)
        
        XSAtime = squeeze(pair.rNB.AStime(rundex(gooddir),:,:));
        XSBtime = squeeze(pair.rNB.BStime(rundex(gooddir),:,:));
        randXSAtime = squeeze(pair.rNB.randAStime(rundex(gooddir),:,:));
        randXSBtime = squeeze(pair.rNB.randBStime(rundex(gooddir),:,:));
                
        XXAtime = squeeze(pair.rNB.AAtime(rundex(gooddir),:,:));
        XXBtime = squeeze(pair.rNB.BBtime(rundex(gooddir),:,:));
        NNtime = squeeze(pair.rNN.NNtime(rundex(gooddir),:,:));
        NNtimecorr = squeeze(pair.rNN.NNtimecorr(rundex(gooddir),:,:));
        
        NNtimecorrsq = squeeze(pair.rNN.NNtimecorr(rundex(gooddir),:,:)).^2;
         
        
        tmp = mean(XSAtime) > thresh.*std(randXSAtime);
        Adex = find(tmp);
        
        tmp = mean(XSBtime) > thresh.*std(randXSBtime);
        Bdex = find(tmp);
        
        dex = intersect(Adex,Bdex);
        if length(find(dex < 250)) > 0
          disp('wtf')
        end
        
        if length(dex) >= minnumtimepts
          % test the standard model
          prod = mean(XSAtime(:,dex)) .* mean(XSBtime(:,dex));
          pairdata.thresh.prod = [pairdata.thresh.prod prod];
          pairdata.thresh.NNcorr = [pairdata.thresh.NNcorr mean(NNtimecorr(:,dex))];
        
          pairdata.thresh.prodsum = [pairdata.thresh.prodsum mean(prod)];
          pairdata.thresh.NNcorrsum = [pairdata.thresh.NNcorrsum mean(mean(NNtimecorr(:,dex)))];
          
          pairdata.thresh.XX = [pairdata.thresh.XX (mean(XXAtime(:,dex))) (mean(XXBtime(:,dex)))];
          
          pairdata.XSthresh = [pairdata.XSthresh mean(XSAtime(:,dex),2)' mean(XSBtime(:,dex),2)'];
          pairdata.XXthresh = [pairdata.XXthresh mean(XXAtime(:,dex),2)' mean(XXBtime(:,dex),2)'];
          pairdata.NNthresh = [pairdata.NNthresh mean(NNtime(:,dex),2)'];
          
          
          pairdata.thresh.XXA = [pairdata.thresh.XXA mean(XXAtime(:,dex),2)'];
          pairdata.thresh.XXB = [pairdata.thresh.XXB mean(XXBtime(:,dex),2)'];
          
          % how did the correlations stack up against random values
          % the 451/250 business is to clip possible values during
          % fixation.
          
          permdex = randperm(451);
          randex = permdex(1:length(dex)) + 250;
          pairdata.thresh.randNNcorr = [pairdata.thresh.randNNcorr mean(NNtimecorrsq(:,randex),2)'];
          
          % and against the peak values?
          [jnk ted] = sort(mean(NNtimecorr),'descend');
          ted = ted(find(ted > 250));
          pairdata.thresh.maxcorr = [pairdata.thresh.maxcorr mean(NNtimecorrsq(:,ted(1:length(dex))),2)'];
          
          % and what was the overlap between their indices?
          compdex = ted(1:length(dex));
          samedex = intersect(compdex,dex);
          pairdata.thresh.frac = [pairdata.thresh.frac length(samedex)/length(dex)];    
          
          % filter dot product (since we're in a loop anyway)
          filtdpm = [mean(mycorr(squeeze(pair.rNB.fA(rundex(gooddir),:,:))',...
            squeeze(pair.rNB.fB(rundex(gooddir),:,:))'))];
          
          filtdps = [std(mycorr(squeeze(pair.rNB.fA(rundex(gooddir),:,:))',...
            squeeze(pair.rNB.fB(rundex(gooddir),:,:))'))];
        
          % correlation of the filters
          pairdata.filtdp.m = [pairdata.filtdp.m filtdpm];
          pairdata.filtdp.s = [pairdata.filtdp.s filtdps];
                  
          % housekeeping
          pairdata.electrodedist = [pairdata.electrodedist pair.electrodedist];
          pairdata.analid = [pairdata.analid pair_ind];
          pairdata.timepts = [pairdata.timepts length(dex)];
        end
        
        
        % and the average trial to trial correlation
        dex = find(ismember(pair.type,prefdirvec(rundex)));
        NNm = [NNm mean(nanmean(pair.rNN.NNcorr(:,find(pair.type == prefdirvec(rundex(gooddir))))))];
        NNs = [NNs std(nanmean(pair.rNN.NNcorr(:,find(pair.type == prefdirvec(rundex(gooddir))))))];
        
      end % by-direction thresholding
       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Other correlational measures %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
      % spike-count correlation
      prs = pair.rSC.prs(rundex,:);
     
      pairdata.rSC.prs.m = [pairdata.rSC.prs.m nanmean(prs,2)'];
      pairdata.rSC.prs.s = [pairdata.rSC.prs.s nanstd(prs,[],2)'];
      
      fix = pair.rSC.fix(rundex,:);
      
      pairdata.rSC.fix.m = [pairdata.rSC.fix.m nanmean(fix,2)'];
      pairdata.rSC.fix.s = [pairdata.rSC.fix.s nanstd(fix,[],2)'];
      
      sig = ttest(pair.rSC.prs(rundex,:),pair.rSC.randprs(rundex,:),.05/50,'right',2);
      pairdata.rSC.prs.sig = [pairdata.rSC.prs.sig sig'];
      
      
      sig = ttest(pair.rSC.fix(rundex,:),pair.rSC.randfix(rundex,:),.05/50,'right',2);
      pairdata.rSC.fix.sig = [pairdata.rSC.fix.sig sig'];

      
      % post-filtered correlation (average, for relevant trials)
            
      pairdata.NNcorr.m = [pairdata.NNcorr.m NNm];
      pairdata.NNcorr.s = [pairdata.NNcorr.s NNs];
      
      
      %%%%%%%%%%%%%%%%
      % housekeeping %
      %%%%%%%%%%%%%%%%
      %pairdata.electrodedist = [pairdata.electrodedist repmat(pair.electrodedist,[1 length(rundex)])];
      
    end % were there significant filters
    
    if ~isempty(nrundex) % were there any non-significant filters?
    
      for gooddir = 1:length(nrundex)
        
        XSAtime = squeeze(pair.rNB.AStime(nrundex(gooddir),:,:));
        XSBtime = squeeze(pair.rNB.BStime(nrundex(gooddir),:,:));
        randXSAtime = squeeze(pair.rNB.randAStime(nrundex(gooddir),:,:));
        randXSBtime = squeeze(pair.rNB.randBStime(nrundex(gooddir),:,:));
                
        XXAtime = squeeze(pair.rNB.AAtime(nrundex(gooddir),:,:));
        XXBtime = squeeze(pair.rNB.BBtime(nrundex(gooddir),:,:));
        NNtime = squeeze(pair.rNN.NNtime(nrundex(gooddir),:,:));
        
        tmp = mean(XSAtime) > thresh.*std(randXSAtime);
        Adex = find(tmp);
        
        tmp = mean(XSBtime) > thresh.*std(randXSBtime);
        Bdex = find(tmp);
        
        dex = intersect(Adex,Bdex);
        if length(find(dex < 250)) > 0
          disp('wtf')
        end
        
        if ~isempty(dex)
          pairdata.uNNthresh = [pairdata.uNNthresh mean(NNtime(:,dex),2)'];
         end
              
      end % by-direction thresholding
      
      
    end % non-significant filters
      
  end % did we run anything at all?

  disp(num2str(pair_ind))
  % clear out the pair's data (for memory reasons)
  clear(sprintf('pair%.3d',pair_ind));

end % aggregation

% fix an error in the fixation rSC
pairdata.rSC.fix.sig(find(isnan(pairdata.rSC.fix.sig))) = 0;


% code to use t-tests instead of thresholded SD
%tmp = ttest(XSAtime,randXSAtime,.05/(50*701),'right');
%tmp = ttest(XSBtime,randXSBtime,.05/(50*701),'right');


ps = [mean(pairdata.thresh.prodsum(1:2)) pairdata.thresh.prodsum(3) ...
  mean(pairdata.thresh.prodsum(4:6)) pairdata.thresh.prodsum(7) ...
  mean(pairdata.thresh.prodsum(8:12)) pairdata.thresh.prodsum(13) ...
  mean(pairdata.thresh.prodsum(14:17)) pairdata.thresh.prodsum(18:22)];

nn = [mean(pairdata.thresh.NNcorrsum(1:2)) pairdata.thresh.NNcorrsum(3) ...
  mean(pairdata.thresh.NNcorrsum(4:6)) pairdata.thresh.NNcorrsum(7) ...
  mean(pairdata.thresh.NNcorrsum(8:12)) pairdata.thresh.NNcorrsum(13) ...
  mean(pairdata.thresh.NNcorrsum(14:17)) pairdata.thresh.NNcorrsum(18:22)];

dp = [mean(pairdata.filtdp.m(1:2)) pairdata.filtdp.m(3) ...
  mean(pairdata.filtdp.m(4:6)) pairdata.filtdp.m(7) ...
  mean(pairdata.filtdp.m(8:12)) pairdata.filtdp.m(13) ...
  mean(pairdata.filtdp.m(14:17)) pairdata.filtdp.m(18:22)];

ed = [mean(pairdata.electrodedist(1:2)) pairdata.electrodedist(3) ...
  mean(pairdata.electrodedist(4:6)) pairdata.electrodedist(7) ...
  mean(pairdata.electrodedist(8:12)) pairdata.electrodedist(13) ...
  mean(pairdata.electrodedist(14:17)) pairdata.electrodedist(18:22)];
