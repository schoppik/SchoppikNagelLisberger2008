%
%
%              by: david schoppik
%            date: 9-28-05
%         purpose: to parse a single direction of trials for analysis
%           usage: parsetrials(experiment,prefdir (a vector is ok),unitnum,
%                  frtype (1 = spikes, 2 = firing rate), binsize
%
%
%
%
%          N.B. will skip trials with no spikes!
%

function [fr,epos,evel,names,type] = parsetrials(experiment,prefdirvec,unitnum,postwin)

if nargin == 0
  help parsetrials
  return
end

%%%%%%%%%%%%%%%%%%%%%%%
%  user's  variables  %
%%%%%%%%%%%%%%%%%%%%%%%

binsize = 1;
frtype = 1;

% how should we align the trials?
% 1. start of pursuit
% 2. start of motion
alignment = 2;

% should we cut saccades - if using the velocity traces?
cutsacs = 1;

% define the filter for the eye movement
[num,denom] = butter(2,25/500);
filter_traces = 1;

% and the smoothing window for the spikes
smfun = gaussian([0 20 1000],-100:100);

% to pull out a particular speed (not all trials)
prefspeed = 20;

% pre/post windows
prewin = 499;
%postwin = 750;

% which trace? 1 = position 2 = velocity 3 = error (target - eye)
%trace = 2;

%%%%%%%%%%%%%%%%%%%%%%%
%  other   variables  %
%%%%%%%%%%%%%%%%%%%%%%%
fr = [];
e = [];

poscal = 40;

%%%%%%%%%%%
%  parse  %
%%%%%%%%%%%

% read in the files
numtrials = 0;
files=dir(pwd);

for i=1:length(files)
  if strcmp(files(i).name(min(length(files(i).name),9)),experiment) && ~strcmp(files(i).name(2:3),'DS')
    numtrials = numtrials +1;
  end
end

ind = 1;
for i=1:numtrials

  % read in the trial data
  trialdata = readcxdata(strcat(getlastdir(pwd),experiment,sprintf('.%04d',i)));

  % if the trial is not the right direction
  if isempty(intersect(str2num(trialdata.trialname(8:10)),prefdirvec));
    continue;
  else
    prefdir = str2num(trialdata.trialname(8:10));
  end

  % if the trial is not the preferred speed
  if length(trialdata.trialname) > 10
    if str2num(trialdata.trialname(12:14)) ~= prefspeed
      continue;
    end
  end

  % now that all trials are marked, there should be no trials w/o a saccade
  % mark except the perfect ones
  if isempty(trialdata.mark1)
    disp(['no mark 1 on trial ' num2str(i)])
    continue;
  end
  
  % if the trial is already defined as a reject
  if trialdata.mark1(1) == -1
    continue;
  end

  % if there were no spikes
  %if trialdata.sortedSpikes{unitnum} < 1
  %  disp('no spikes')
  %  continue;
  %end

  % if there were no saccades marked
  if length(trialdata.mark1) > 1 & isempty(trialdata.mark2)
    continue;
  end

  % or if marking was otherwise screwed up
  if length(trialdata.mark1) < length(trialdata.mark2)
    continue;
  end

  % align the trials
  switch alignment
    case 1
      start = trialdata.mark1(max(find(trialdata.mark1)));
    case 2
      start = floor(trialdata.other(7));
    case 3;
      start = floor(trialdata.other(6)) - (prewin-200-1);
  end

  % if the unit wasn't being recorded
  if trialdata.sortedSpikes{unitnum} == -99
    continue;
  end

  spvec = zeros(1,length(trialdata.data));

  % if the unit responded at all
  if trialdata.sortedSpikes{unitnum}(1) > 0
    spvec(ceil(trialdata.sortedSpikes{unitnum})) = 1;
  end
    
  switch frtype
    case 1
      % recreate the binary spike train
      fr(ind,:) = binspikedata(spvec(start-prewin:start+postwin),binsize);
    case 2
      frate = conv(spvec(start-prewin:min(postwin+start,length(trialdata.data))),smfun);
      fr(ind,:) = bindata(frate(101:101+prewin+postwin),binsize);
  end
  
  % get the eye movement
  hepos = filtfilt(num,denom,trialdata.data(1,:)./poscal);
  vepos = filtfilt(num,denom,trialdata.data(2,:)./poscal);
  repos = sqrt((hepos.^2 + vepos.^2));
  
  heye = filtfilt(num,denom,digitaldiff(trialdata.data(1,:)./poscal));
  veye = filtfilt(num,denom,digitaldiff(trialdata.data(2,:)./poscal));

  if cutsacs
    if max(trialdata.mark1) > length(trialdata.data)
      disp(num2str(i))
      continue;
    else
      sacstart = [trialdata.mark1];
      sacend = [trialdata.mark2];
      
      hevel = removesacs(heye,sacstart,sacend);
      vevel = removesacs(veye,sacstart,sacend);
      revel = sqrt((hevel.^2 + vevel.^2));
    end
  end


  epos(ind,:,1) =  bindata(hepos(start-prewin:start+postwin),binsize);
  epos(ind,:,2) =  bindata(vepos(start-prewin:start+postwin),binsize);
  evel(ind,:,1) = bindata(hevel(start-prewin:start+postwin),binsize);
  evel(ind,:,2) = bindata(vevel(start-prewin:start+postwin),binsize);

  %   switch prefdir
%     case {360,180}
%       epos(ind,:) = bindata(hepos(start-prewin:start+postwin),binsize);
%       evel(ind,:) = bindata(hevel(start-prewin:start+postwin),binsize);
%     case {90 270}
%       epos(ind,:) = bindata(vepos(start-prewin:start+postwin),binsize);
%       evel(ind,:) = bindata(vevel(start-prewin:start+postwin),binsize);
%     otherwise
%       epos(ind,:) = bindata(repos(start-prewin:start+postwin),binsize);
%       evel(ind,:) = bindata(revel(start-prewin:start+postwin),binsize);
%   end

  type(ind) = prefdir;
  names(ind,:) = strcat(getlastdir(pwd),experiment,sprintf('.%04d',i));
  ind = ind+1;
end

function [eye] = removesacs(eye,sacstart,sacend);

offset = 5;

for i = 1:length(sacend)
  % interpolate between the start of the saccade and the end
  ind_start = max(1,sacstart(i)-offset);
  ind_end = min(sacend(i)+offset,length(eye));
  
  if sacstart(i)-offset < 1
    ind_start = 1;
    val_start = 0;
  else
    ind_start = sacstart(i)-offset;
    val_start = eye(ind_start);
  end
  
  if sacend(i)+offset > length(eye)
    ind_end = length(eye);
    val_end = val_start;
  else
    ind_end = sacend(i) + offset;
    val_end = eye(ind_end);
  end
  
  if val_start == val_end
    eye(ind_start:ind_end) = val_start*ones(1,ind_end-ind_start+1);

  else
    stepsize = (val_end-val_start)/(ind_end-ind_start);
    eye(ind_start:ind_end) = val_start:stepsize:val_end;
  end

end
