%
%
%
%         by: david schoppik
%       date: 4/16/07
%    purpose: to measure the relationship between a pair of neurons
%             in this iteration, we'll filter one spike trace, so it is
%             reasonably close to the "eye movement" filtering that we've
%             performed.  we'll keep the spike-count correlation as well so
%             as to have a baseline to compare it to.

function pair = nbnv_final(pair);
prefdirvec = 45:45:360;
minSpikesPerTrial = 4;
numboot = 50;
use_fraction = .4;

shortrun = 0;

if ~shortrun
  % remove old shit
  pair = rmfield(pair,'rNN');
  % filters
  pair.rNN.fN = single(zeros(8,numboot,701));
  pair.rNN.randfN = single(zeros(8,numboot,701));

  % predictions
  pair.rNN.NN = single(nan(numboot,length(pair.type)));
  pair.rNN.NNtime = single(zeros(8,numboot,701));

  % normalized predictions
  pair.rNN.NNcorr = single(nan(numboot,length(pair.type)));
  pair.rNN.NNtimecorr = single(zeros(8,numboot,701));
  
  pair = rmfield(pair,'rNB');
  % filters
  pair.rNB.fA = single(zeros(8,numboot,701));
  pair.rNB.fB = single(zeros(8,numboot,701));
  pair.rNB.randfA = single(zeros(8,numboot,701));
  pair.rNB.randfB = single(zeros(8,numboot,701));

  % correlation and covariation across trials
  pair.rNB.AS = single(nan(numboot,length(pair.type)));
  pair.rNB.BS = single(nan(numboot,length(pair.type)));
  pair.rNB.AStime = single((zeros(8,numboot,701)));
  pair.rNB.BStime = single((zeros(8,numboot,701)));
  pair.rNB.AA = single(nan(numboot,length(pair.type)));
  pair.rNB.BB = single(nan(numboot,length(pair.type)));
  pair.rNB.AAtime = single((zeros(8,numboot,701)));
  pair.rNB.BBtime = single((zeros(8,numboot,701)));

  % and for the examplepair
  pair.rNB.randAStime = single((zeros(8,numboot,701)));  
  pair.rNB.randBStime = single((zeros(8,numboot,701)));
  
  % spike count correlation
  pair = rmfield(pair,'rSC');
  pair.rSC.fix = single(zeros(8,numboot));
  pair.rSC.prs = single(zeros(1,8));
  pair.rSC.randfix = single(zeros(8,numboot));
  pair.rSC.randprs = single(zeros(8,numboot));
end

if shortrun
  % and when we don't want to do all that bootstrapping stuff
  pair = rmfield(pair,'fast');
  pair.fast.AAcov = zeros(8,701);
  pair.fast.BBcov = zeros(8,701);
  pair.fast.AScorr = zeros(8,701);
  pair.fast.BScorr = zeros(8,701);
  pair.fast.cov = zeros(8,701);
  
  pair = rmfield(pair,'short');
  pair.short.AA = zeros(1,8);
  pair.short.BB = zeros(1,8);
  pair.short.AS = zeros(1,8);
  pair.short.BS = zeros(1,8);
  pair.short.NN = zeros(1,8);
  pair.short.rSC = zeros(1,8);
end


for dir = 1:length(prefdirvec)
  dex = find(pair.type == prefdirvec(dir));
  numtrials = length(dex);

  if sum(sum(pair.spikes(dex,:,1)))/length(dex) >= minSpikesPerTrial &...
      sum(sum(pair.spikes(dex,:,2)))/length(dex) >= minSpikesPerTrial
    pair.run(dir) = 1;
    bootkeep = round(length(dex)*use_fraction);


    spA = pair.spikes(dex,300:end,1);
    spB = pair.spikes(dex,300:end,2);
    E = abs(pair.evel(dex,300:end));

    % get the variability
    spAvar = spA - repmat(mean(spA),[size(spA,1) 1]);
    spBvar = spB - repmat(mean(spB),[size(spA,1) 1]);
    Evar = E - repmat(mean(E),[size(spA,1) 1]);

    if ~shortrun
      % since we're bootstrapping, we first fill arrays first with randomly
      % shuffled trials, so that we can use vector operations later
      fulle3D = single(zeros([fliplr(size(E)) numboot]));
      fullAspikes3D = single(zeros([fliplr(size(spA)) numboot]));
      fullBspikes3D = single(zeros([fliplr(size(spB)) numboot]));

      % permute
      [jnk perm] = sort(rand(numboot,numtrials),2);

      % fill the arrays
      % I couldn't find an easy non-loop way to do this, though it no
      % doubt exists with some tricky reshaping.
      for draw = 1:numboot
        fulle3D(:,:,draw) = E(perm(draw,:),:)';
        fullAspikes3D(:,:,draw) = spA(perm(draw,:),:)';
        fullBspikes3D(:,:,draw) = spB(perm(draw,:),:)';
      end

      % subtract the mean in time -- we're doing the analysis on the covariance
      vare3D = fulle3D - repmat(mean(fulle3D,2),[1,numtrials,1]);
      varAspikes3D = fullAspikes3D - repmat(mean(fullAspikes3D,2),[1,numtrials,1]);
      varBspikes3D = fullBspikes3D - repmat(mean(fullBspikes3D,2),[1,numtrials,1]);

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
      varSA_filt = varAspikes3D(:,useindex,:);
      varSB_filt = varBspikes3D(:,useindex,:);

      % matrices of variability for testing
      varE_test = vare3D(:,testindex,:);
      varSA_test = varAspikes3D(:,testindex,:);
      varSB_test = varBspikes3D(:,testindex,:);

      % randomly permute the matrix of variability for filter generation to
      % test the statistical significance of the filter

      varSA_rand = single(zeros(size(varSA_filt)));
      varSB_rand = single(zeros(size(varSB_filt)));

      [jnk filtperm] = sort(rand(numboot,numfilttrials),2);
      for draw = 1:numboot
        varSA_rand(:,:,draw) = varSA_filt(:,filtperm(draw,:),draw);
        varSB_rand(:,:,draw) = varSB_filt(:,filtperm(draw,:),draw);
      end

      %%%%%%%%%%%
      % filters %
      %%%%%%%%%%%

      % generate the filters
      varAsta = mean(getfilter(varSA_filt,varE_filt),2);
      varBsta = mean(getfilter(varSB_filt,varE_filt),2);

      shufvarAsta = mean(getfilter(varSA_rand,varE_filt),2);
      shufvarBsta = mean(getfilter(varSB_rand,varE_filt),2);

      % generate the spike-spike filter
      varNsta = mean(getfilter(varSA_filt,varSB_filt),2);
      shufvarNsta = mean(getfilter(varSA_filt,varSB_rand),2);

      % we'll have to interpolate the zeroth value here, if appropriate
      if ~pair.electrodedist
        for draw = 1:numboot
          varNsta(:,1,draw) = interp_zero_lag(varNsta(:,1,draw));
        end
      end

      % save the filters
      pair.rNB.fA(dir,:,:) = squeeze(single(varAsta))';
      pair.rNB.fB(dir,:,:) = squeeze(single(varBsta))';
      pair.rNN.fN(dir,:,:) = squeeze(single(varNsta))';

      pair.rNB.randfA(dir,:,:) = squeeze(single(shufvarAsta))';
      pair.rNB.randfB(dir,:,:) = squeeze(single(shufvarBsta))';
      pair.rNN.randfN(dir,:,:) = squeeze(single(shufvarNsta))';
      %%%%%%%%%%%%%%%%%%%%
      % test the filters %
      %%%%%%%%%%%%%%%%%%%%

      % use the filters to make predictions about the variability on each trial
      varEA_pred = myconv(varSA_test,repmat(varAsta,[1,size(varSA_test,2),1]));
      varEB_pred = myconv(varSB_test,repmat(varBsta,[1,size(varSB_test,2),1]));

      shufvarEA_pred = myconv(varSA_test,repmat(shufvarAsta,[1,size(varSA_test,2),1]));
      shufvarEB_pred = myconv(varSB_test,repmat(shufvarBsta,[1,size(varSB_test,2),1]));


      % test each filter against the shuffled data
      for draw = 1:numboot

        %%%%%%
        % NB %
        %%%%%%
        predicted_Avar = varEA_pred(:,:,draw); % predicted
        predicted_Bvar = varEB_pred(:,:,draw); % predicted
        actual_var = varE_test(:,:,draw); % actual
        
        shufpredicted_Avar = shufvarEA_pred(:,:,draw); % predicted
        shufpredicted_Bvar = shufvarEB_pred(:,:,draw); % predicted
        
        % keep track of how we do for each trial
        pair.rNB.AS(draw,dex(perm(draw,testindex))) = single(mycorr(predicted_Avar,actual_var));
        pair.rNB.BS(draw,dex(perm(draw,testindex))) = single(mycorr(predicted_Bvar,actual_var));

        pair.rNB.AStime(dir,draw,:) = single(mycorr(predicted_Avar',actual_var'));
        pair.rNB.BStime(dir,draw,:) = single(mycorr(predicted_Bvar',actual_var'));
      
        pair.rNB.AA(draw,dex(perm(draw,testindex))) = single(mycov(predicted_Avar,predicted_Avar));
        pair.rNB.BB(draw,dex(perm(draw,testindex))) = single(mycov(predicted_Bvar,predicted_Bvar));

        pair.rNB.AAtime(dir,draw,:) = single(mycov(predicted_Avar',predicted_Avar'));
        pair.rNB.BBtime(dir,draw,:) = single(mycov(predicted_Bvar',predicted_Bvar'));

        pair.rNB.randAStime(dir,draw,:) = single(mycorr(shufpredicted_Avar',actual_var'));
        pair.rNB.randBStime(dir,draw,:) = single(mycorr(shufpredicted_Bvar',actual_var'));
 
        
        %%%%%%
        % NN %
        %%%%%%

        pair.rNN.NN(draw,dex(perm(draw,testindex))) = single(mycov(predicted_Avar,predicted_Bvar));
        pair.rNN.NNtime(dir,draw,:) = single(mycov(predicted_Avar',predicted_Bvar'));

        pair.rNN.NNcorr(draw,dex(perm(draw,testindex))) = single(mycorr(predicted_Avar,predicted_Bvar));
        pair.rNN.NNtimecorr(dir,draw,:) = single(mycorr(predicted_Avar',predicted_Bvar'));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spike count correlation %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        bootdex = perm(draw,1:bootkeep);
        shufdex = bootdex(randperm(length(bootdex)));

        % standard noise correlation across trials, split by fixation and
        % behavior
        pair.rSC.fix(dir,draw) = mycorr(squeeze(sum(pair.spikes(dex(bootdex),1:500,:),2)));
        pair.rSC.prs(dir,draw) = mycorr(squeeze(sum(pair.spikes(dex(bootdex),501:end,:),2)));
        pair.rSC.randfix(dir,draw) = mycorr([sum(pair.spikes(dex(bootdex),1:500,1),2),sum(pair.spikes(dex(shufdex),1:500,2),2)]);
        pair.rSC.randprs(dir,draw) = mycorr([sum(pair.spikes(dex(bootdex),501:end,1),2),sum(pair.spikes(dex(shufdex),501:end,2),2)]);

         
      end

    end
    
    if shortrun
      % without all that fancy-pants bootstrapping action
      fA = mean(getfilter(spAvar',Evar'),2);
      fB = mean(getfilter(spBvar',Evar'),2);

      EhatA = myconv(spAvar',repmat(fA,[1 size(spA,1)]));
      EhatB = myconv(spBvar',repmat(fB,[1 size(spA,1)]));
      
      pair.fast.AAcov(dir,:) = (mycov(EhatA',EhatA'));
      pair.fast.BBcov(dir,:) = (mycov(EhatB',EhatB'));
      
      pair.fast.AScorr(dir,:) = (mycorr(EhatA',Evar));
      pair.fast.BScorr(dir,:) = (mycorr(EhatB',Evar));
      
      pair.fast.cov(dir,:) = (mycov(EhatA',EhatB'));
      
      pair.short.AA(dir) = fitskewnorm(mycov(EhatA,EhatA),1);
      pair.short.BB(dir) = fitskewnorm(mycov(EhatB,EhatB),1);
      pair.short.AS(dir) = mean(mycorr(EhatA,Evar'));
      pair.short.BS(dir) = mean(mycorr(EhatB,Evar'));
      pair.short.NN(dir) = fitskewnorm(mycov(EhatA,EhatB),1);
      pair.short.rSC(dir) = mycorr(squeeze(sum(pair.spikes(dex,501:end,:),2)));
    end
    
  else
    pair.run(dir) = 0;

  end
end


%  pair.rNN.NNcov = single(zeros(8,numboot));
%  pair.rNN.NNtimecov = single(zeros(8,numboot,701));


%         pair.rNB.AS(dir,draw) = median(mycorr(predicted_Avar,actual_var));
%         pair.rNB.BS(dir,draw) = median(mycorr(predicted_Bvar,actual_var));
% 
%         pair.rNB.Atime(dir,draw,:) = single(mycorr(predicted_Avar',actual_var'));
%         pair.rNB.Btime(dir,draw,:) = single(mycorr(predicted_Bvar',actual_var'));
% 
%         pair.rNB.AA(dir,draw) = median(mycov(predicted_Avar,predicted_Avar));
%         pair.rNB.BB(dir,draw) = median(mycov(predicted_Bvar,predicted_Bvar));
% 
%         pair.rNB.Atimecov(dir,draw,:) = single(mycov(predicted_Avar',predicted_Avar'));
%         pair.rNB.Btimecov(dir,draw,:) = single(mycov(predicted_Bvar',predicted_Bvar'));

%        pair.rNN.NNcov(dir,draw) = median(mycov(predicted_Avar,predicted_Bvar));
%        pair.rNN.NNtimecov(dir,draw,:) = single(mycov(predicted_Avar',predicted_Bvar'));


%pair.rNB.randA = single(zeros(8,numboot));
%pair.rNB.randB = single(zeros(8,numboot));
%pair.rNB.randAtime = single(zeros(8,numboot,701));
%pair.rNB.randBtime = single(zeros(8,numboot,701));


%   shufvarEA_pred = myconv(varSA_test,repmat(shufvarAsta,[1,size(varSA_test,2),1]));
%   shufvarEB_pred = myconv(varSB_test,repmat(shufvarBsta,[1,size(varSB_test,2),1]));

% and now the N-N
%   varNN_pred = myconv(varSA_test,repmat(varNsta,[1,size(varSB_test,2),1]));
%   shufvarNN_pred = myconv(varSA_test,repmat(shufvarNsta,[1,size(varSB_test,2),1]));

%     shufpredicted_Avar = shufvarEA_pred(:,:,draw); % predicted
%     shufpredicted_Bvar = shufvarEB_pred(:,:,draw); % predicted
%      pair.rNB.randA(dir,draw) = median(mycorr(shufpredicted_Avar,actual_var));
%      pair.rNB.randB(dir,draw) = single(corr2(shufpredicted_Bvar(:),actual_var(:)));
%      pair.rNB.randAtime(dir,draw,:) = single(mycorr(shufpredicted_Avar',actual_var'));
%      pair.rNB.randBtime(dir,draw,:) = single(mycorr(shufpredicted_Bvar',actual_var'));


%       %%%%%%
%       % NN %
%       %%%%%%
%
%       predicted_Nvar = varNN_pred(:,:,draw); % predicted
%       shufpredicted_Nvar = shufvarNN_pred(:,:,draw);
%       actual_Nvar = varSB_test(:,:,draw); % actual
%
%       pair.rNN.corr(dir,draw) = single(corr2(predicted_Nvar(:),actual_Nvar(:)));
%       pair.rNN.randcorr(dir,draw) = single(corr2(shufpredicted_Nvar(:),actual_Nvar(:)));
%       pair.rNN.time(dir,draw,:) = single(mycorr(predicted_Nvar',actual_Nvar'));
%       pair.rNN.randtime(dir,draw,:) = single(mycorr(shufpredicted_Nvar',actual_Nvar'));
%
%       pair.rNN.conv(dir,draw,:) = single(mean(myconv(predicted_Avar,predicted_Bvar),2));
%


%     % filter one of the spike trains, so that it has the same frequency
%     % content as the eye movement
%     spBsmooth = zeros(size(spB));
%
%     for i = 1:length(dex)
%       spBsmooth(i,:) = filtfilt(num,denom,spB(i,:));
%     end
%     % prepare for the 3d run
%




% get the coherence
%pair.rNN.coh(dir,:) = mycohere(varSA_filt,varSB_filt);
%pair.rNN.paircoh(dir,:) = mycohere(varEA_pred,varEB_pred);

%pair.rNN.randcoh(dir,:) = mycohere(varSA_rand,varSB_filt);
%pair.rNN.randpaircoh(dir,:) = mycohere(shufvarEA_pred,varEB_pred);


% compare the post-filtered prediction
%pair.rNN.paircorr(dir,draw) = mean(single(mycorr(predicted_Avar,predicted_Bvar)));
%pair.rNN.randpaircorr(dir,draw) = mean(single(mycorr(shufpredicted_Avar,shufpredicted_Bvar)));
%pair.rNN.pairtime(dir,draw,:) = single(mycorr(predicted_Avar',predicted_Bvar'));
%pair.rNN.randpairtime(dir,draw,:) = single(mycorr(shufpredicted_Avar',shufpredicted_Bvar'));


% pair.rNN.corr = single(zeros(8,numboot));
% pair.rNN.randcorr = single(zeros(8,numboot));
% pair.rNN.time = single(zeros(8,numboot,701));
% pair.rNN.randtime = single(zeros(8,numboot,701));
% pair.rNN.coh = single(zeros(8,numboot));
% pair.rNN.paircoh = single(zeros(8,numboot));
% pair.rNN.randcoh = single(zeros(8,numboot));
% pair.rNN.randpaircoh = single(zeros(8,numboot));
% pair.rNN.paircorr = single(zeros(8,numboot));
% pair.rNN.randpaircorr = single(zeros(8,numboot));
% pair.rNN.pairtime = single(zeros(8,numboot,701));
% pair.rNN.randpairtime = single(zeros(8,numboot,701));

% pair.rNB.fonlyA = single(zeros(8,701));
% pair.rNB.fonlyB = single(zeros(8,701));
% pair.rNB.fboth = single(zeros(8,701));
%
% pair.rNB.sep = single(zeros(8,numboot));
% pair.rNB.tog = single(zeros(8,numboot));
% pair.rNB.septime = single(zeros(8,numboot,701));
% pair.rNB.togtime = single(zeros(8,numboot,701));
%
% pair.rNB.fonlyA = single(zeros(8,701));
% pair.rNB.fonlyB = single(zeros(8,701));
% pair.rNB.fboth = single(zeros(8,701));
%
% pair.rNB.fsum = single(zeros(8,701));
% pair.rNB.sum = single(zeros(8,numboot));
% pair.rNB.sumtime = single(zeros(8,numboot,701));

% now, for the Dan analysis
% spAc = countspikedata(spA,window_size);
% spBc = countspikedata(spB,window_size);

% sponlyA = zeros(size(spA));
% sponlyB = zeros(size(spA));
% sponlyP = zeros(size(spA));

% find the times A fired and B did not
% sponlyA(find(spAc & ~spBc)) = spAc(find(spAc & ~spBc));
% sponlyB(find(spBc & ~spAc)) = spBc(find(spBc & ~spAc));

% find the pairs
% spboth = spAc + spBc;
% spboth(find(spboth == 1)) = 0;
%     fullonlyAspikes3D = single(zeros([fliplr(size(sponlyA)) numboot]));
%     fullonlyBspikes3D = single(zeros([fliplr(size(sponlyB)) numboot]));
%     fullbothspikes3D = single(zeros([fliplr(size(spboth)) numboot]));
%     fullsumspikes3D = single(zeros([fliplr(size(spA + spB)) numboot]));
%       fullonlyAspikes3D(:,:,draw) = sponlyA(perm(draw,:),:)';
%       fullonlyBspikes3D(:,:,draw) = sponlyB(perm(draw,:),:)';
%       fullbothspikes3D(:,:,draw) = spboth(perm(draw,:),:)';
%     varonlyAspikes3D = fullonlyAspikes3D - repmat(mean(fullonlyAspikes3D,2),[1,numtrials,1]);
%     varonlyBspikes3D = fullonlyBspikes3D - repmat(mean(fullonlyBspikes3D,2),[1,numtrials,1]);
%     varbothspikes3D = fullbothspikes3D - repmat(mean(fullbothspikes3D,2),[1,numtrials,1]);
%     varsumspikes3D = varAspikes3D + varBspikes3D;



%     varSonlyA_test = varonlyAspikes3D(:,testindex,:);
%     varSonlyB_test = varonlyBspikes3D(:,testindex,:);
%     varSboth_test = varbothspikes3D(:,testindex,:);
%     varSsum_test = varsumspikes3D(:,testindex,:);
%

%     varSonlyA_filt = varonlyAspikes3D(:,useindex,:);
%     varSonlyB_filt = varonlyBspikes3D(:,useindex,:);
%     varSboth_filt = varbothspikes3D(:,useindex,:);
%     varSsum_filt = varsumspikes3D(:,useindex,:);


%     varonlyAsta = mean(myconv(varSonlyA_filt,varE_filt),2);
%     varonlyBsta = mean(myconv(varSonlyB_filt,varE_filt),2);
%     varbothsta = mean(myconv(varSboth_filt,varE_filt),2);
%     varsumsta = mean(myconv(varSsum_filt,varE_filt),2);


%     pair.rNB.fonlyA(dir,:) = single(mean(varonlyAsta,3));
%     pair.rNB.fonlyB(dir,:) = single(mean(varonlyBsta,3));
%     pair.rNB.fboth(dir,:) = single(mean(varbothsta,3));
%     pair.rNB.fsum(dir,:) = single(mean(varsumsta,3));

%     varEonlyA_pred = myconv(varSonlyA_test,repmat(varonlyAsta,[1,size(varSonlyA_test,2),1]));
%     varEonlyB_pred = myconv(varSonlyB_test,repmat(varonlyBsta,[1,size(varSonlyB_test,2),1]));
%     varEboth_pred = myconv(varSboth_test,repmat(varbothsta,[1,size(varSboth_test,2),1]));
%     varEsum_pred = myconv(varSsum_test,repmat(varsumsta,[1,size(varSsum_test,2),1]));

% now sum the predictions to test things against one another
%     varEsep = varEA_pred + varEB_pred;
%     varEtog = varEonlyA_pred + varEonlyB_pred + varEboth_pred;
%

%       predicted_sep = varEsep(:,:,draw);
%       predicted_tog = varEtog(:,:,draw);
%       predicted_sum = varEsum_pred(:,:,draw);
%

% for the separate/together analysis
%       pair.rNB.sep(dir,draw) = single(corr2(predicted_sep(:),actual_var(:)));
%       pair.rNB.tog(dir,draw) = single(corr2(predicted_tog(:),actual_var(:)));
%
%       pair.rNB.septime(dir,draw,:) = single(mycorr(predicted_sep',actual_var'));
%       pair.rNB.togtime(dir,draw,:) = single(mycorr(predicted_tog',actual_var'));
%
%       % for the sum
%       pair.rNB.sum(dir,draw) = single(mycorr(predicted_sum(:),actual_var(:)));
%       pair.rNB.sumtime(dir,draw,:) = single(mycorr(predicted_sum',actual_var'));
%
% for the NN
%       predicted_Nvar = countspikedata(varNN_pred(:,:,draw)',window_size)';
%       shufpredicted_Nvar = countspikedata(shufvarNN_pred(:,:,draw)',window_size)';
%       actual_Nvar = countspikedata(varSB_test(:,:,draw)',window_size)';
%
%       pair.rNN.corr(dir,draw) = single(corr2(predicted_Nvar(:),actual_Nvar(:)));
%       pair.rNN.randcorr(dir,draw) = single(corr2(shufpredicted_Nvar(:),actual_Nvar(:)));
%
%       pair.rNN.time(dir,draw,:) = single(mycorr(predicted_Nvar',actual_Nvar'));
%       pair.rNN.randtime(dir,draw,:) = single(mycorr(shufpredicted_Nvar',actual_Nvar'));




%fullsumspikes3D(:,:,draw) = spsum(perm(draw,:),:)';

%     for boot = 1:numboot
%
%       % randomly draw from the trials
%       bootdex = randperm(length(dex));
%
%       % only use a fraction of them
%       testdex = bootdex(bootkeep+1:end);
%       bootdex = bootdex(1:bootkeep);
%
%       % shuffle the data
%       shufdex = bootdex(randperm(length(bootdex)));
%       shufeyedex = bootdex(randperm(length(bootdex)));
%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % Spike count correlation %
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % standard noise correlation across trials, split by fixation and
%       % behavior
%       pair.rSC.fix(dir,boot) = mycorr(squeeze(sum(pair.spikes(dex(bootdex),1:500,:),2)));
%       pair.rSC.prs(dir,boot) = mycorr(squeeze(sum(pair.spikes(dex(bootdex),501:end,:),2)));
%       pair.rSC.randfix(dir,boot) = mycorr([sum(pair.spikes(dex(bootdex),1:500,1),2),sum(pair.spikes(dex(shufdex),1:500,2),2)]);
%       pair.rSC.randprs(dir,boot) = mycorr([sum(pair.spikes(dex(bootdex),501:end,1),2),sum(pair.spikes(dex(shufdex),501:end,2),2)]);
%
%       % now for appropriate comparison, using the -200:500 of the trial
%       pair.rNN.corr(dir,boot) = mycorr(squeeze(sum(pair.spikes(dex(bootdex),300:end,:),2)));
%       pair.rNN.randcorr(dir,boot) =mycorr([sum(pair.spikes(dex(bootdex),300:end,1),2),sum(pair.spikes(dex(shufdex),300:end,2),2)]);
%
%       % now, test the spikes in time
%       fN = mean(myconv(varspA(bootdex,300:end)',varspB(bootdex,300:end)'),2);
%       shufN = mean(myconv(varspA(bootdex,300:end)',varspB(shufdex,300:end)'),2);
%
%       actsp = varspA(testdex,300:end)';
%       predsp = myconv(varspB(testdex,300:end)',repmat(fN,[1,length(testdex)]));
%       predshufsp = myconv(varspB(testdex,300:end)',repmat(shufN,[1,length(testdex)]));
%
%       pair.rNN.time(dir,boot,:) = mycorr(actsp',predsp');
%       pair.rNN.randtime(dir,boot,:) = mycorr(actsp',predshufsp');
%
%       % get the filters
%       fA = getfilter(varspA(bootdex,300:end)',varE(bootdex,300:end)');
%       fB = getfilter(varspB(bootdex,300:end)',varE(bootdex,300:end)');
%
%       shufA = getfilter(varspA(bootdex,300:end)',varE(shufdex,300:end)');
%       shufB = getfilter(varspB(bootdex,300:end)',varE(shufdex,300:end)');
%
%       % test the filters
%       actE = varE(testdex,300:end)';
%       predEA = myconv(varspA(testdex,300:end)',repmat(fA,[1,length(testdex)]));
%       predEB = myconv(varspB(testdex,300:end)',repmat(fB,[1,length(testdex)]));
%
%       predshufEA = myconv(varspA(testdex,300:end)',repmat(shufA,[1,length(testdex)]));
%       predshufEB = myconv(varspB(testdex,300:end)',repmat(shufB,[1,length(testdex)]));
%
%       pair.rNB.A(dir,boot) = corr2(actE(:),predEA(:));
%       pair.rNB.B(dir,boot) = corr2(actE(:),predEB(:));
%       pair.rNB.randA(dir,boot) = corr2(actE(:),predshufEA(:));
%       pair.rNB.randB(dir,boot) = corr2(actE(:),predshufEB(:));
%
%       % in time
%       pair.rNB.Atime(dir,boot,:) = mycorr(actE',predEA');
%       pair.rNB.Btime(dir,boot,:) = mycorr(actE',predEB');
%       pair.rNB.randAtime(dir,boot,:) = mycorr(actE',predshufEA');
%       pair.rNB.randBtime(dir,boot,:) = mycorr(actE',predshufEA');
%
%     end



%varspA = spA - repmat(mean(spA),[size(spA,1),1]);
%varspB = spB - repmat(mean(spB),[size(spB,1),1]);
%varE = E - repmat(mean(E),[size(E,1),1]);



% now, some summary statistics

% if ~isempty(find(pair.run))
%   fixttest = ttest(pair.rSC.fix(find(pair.run),:),pair.rSC.randfix(find(pair.run),:),.05/100,'both',2);
%   prsttest = ttest(pair.rSC.prs(find(pair.run),:),pair.rSC.randprs(find(pair.run),:),.05/100,'both',2);
%
%   rundex = find(pair.run);
%   evalfix = rundex(find(fixttest));
%   evalprs = rundex(find(prsttest));
%
%   if ~isempty(evalfix)
%     pair.summary.rSC.fix.mean = mean(mean(pair.rSC.fix(find(pair.run),:),2));
%     pair.summary.rSC.fix.std = std(mean(pair.rSC.fix(find(pair.run),:),2));
%     pair.summary.rSC.randfix.mean = mean(mean(pair.rSC.randfix(find(pair.run),:),2));
%     pair.summary.rSC.randfix.std = std(mean(pair.rSC.randfix(find(pair.run),:),2));
%   end
%
%   if ~isempty(evalprs)
%     pair.summary.rSC.prs.mean = mean(mean(pair.rSC.prs(find(pair.run),:),2));
%     pair.summary.rSC.prs.std = std(mean(pair.rSC.prs(find(pair.run),:),2));
%     pair.summary.rSC.randprs.mean = mean(mean(pair.rSC.randprs(find(pair.run),:),2));
%     pair.summary.rSC.randprs.std = std(mean(pair.rSC.randprs(find(pair.run),:),2));
%   end
%
%   pair.summary.rSC.fixttest = evalfix;
%   pair.summary.rSC.prsttest = evalprs;
% else
%   pair.summary.rSC.fixtest = nan; % there was nothing to test
%   pair.summary.rSC.prstest = nan;
% end
%
