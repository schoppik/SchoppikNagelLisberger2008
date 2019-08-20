%
%
%
%
%          by: david schoppik
%        date: 08-14-2006
%     purpose: To go from data sorted by date of experiment to data sorted
%              by recorded unit -- creates "unit???.mat" files
%
%       calls: different variants of the "tbtv" function, which is the main
%              analysis code.  it was easier historically to adapt the tbtv
%              code to run on different data sets than to write one really
%              big function that did everything.
%              parsetrials.m is the main code to read in an experiment worth
%              of data
%              minimize_variance.m will align the data so as to minimize
%              the variability in the eye movements.

tic

% variables
prefdirvec = 45:45:360;
aboveflankvec = circshift(45:45:360,[0 -1]);
belowflankvec = circshift(45:45:360,[0 1]);
data_dir = '/home/schoppik/exp/alldata/';
image_dir = '/home/schoppik/exp/images/';
output_dir = '/home/schoppik/exp/long/';

% just load the pursuit database
load pursdb
uniqueids = intersect(vertcat(p.id),1:601);
allids = vertcat(p.id);

% fixes
p(18).fullid = 'cb040707.d.sig002a';
p(18).localid = 2;
p(115).fullid = 'cb050808.b.sig004a';
p(115).localid = 4;

% for each unique ID
for id = 1:length(uniqueids)
   
  % see the code at the bottom to adapt for trials of different lengths
  % this length is a compromise
  postwin = 500;
  
  % find out how many entries in p there are that match each unit
  matching_files = find(allids == uniqueids(id));
  
  % concatenate the data, if need be -- there are occasionally days where
  % the same unit persisted across different data files/experiments
  switch length(matching_files)
    case 1
      cd(strcat(data_dir,p(matching_files(1)).fullid(1:8)));
      experiment = p(matching_files(1)).fullid(10);
      unitnum = p(matching_files(1)).localid;
      [unit.allspikes,unit.allepos,unit.allevel,unit.allnames,unit.alltype] = parsetrials(experiment,prefdirvec,unitnum,postwin);
    case 2
      cd(strcat(data_dir,p(matching_files(1)).fullid(1:8)));
      experiment = p(matching_files(1)).fullid(10);
      unitnum = p(matching_files(1)).localid;
      [spikes1,epos1,evel1,names1,type1] = parsetrials(experiment,prefdirvec,unitnum,postwin);

      cd(strcat(data_dir,p(matching_files(2)).fullid(1:8)));
      experiment = p(matching_files(2)).fullid(10);
      unitnum = p(matching_files(2)).localid;
      [spikes2,epos2,evel2,names2,type2] = parsetrials(experiment,prefdirvec,unitnum,postwin);
      
      unit.allspikes = [spikes1;spikes2];
      unit.allepos = [epos1;epos2];
      unit.allevel = [evel1;evel2];
      unit.allnames = [names1;names2];
      unit.alltype = [type1 type2];
      
    case 3
      cd(strcat(data_dir,p(matching_files(1)).fullid(1:8)));
      experiment = p(matching_files(1)).fullid(10);
      unitnum = p(matching_files(1)).localid;
      [spikes1,epos1,evel1,names1,type1] = parsetrials(experiment,prefdirvec,unitnum,postwin);

      cd(strcat(data_dir,p(matching_files(2)).fullid(1:8)));
      experiment = p(matching_files(2)).fullid(10);
      unitnum = p(matching_files(2)).localid;
      [spikes2,epos2,evel2,names2,type2] = parsetrials(experiment,prefdirvec,unitnum,postwin);
  
      cd(strcat(data_dir,p(matching_files(3)).fullid(1:8)));
      experiment = p(matching_files(3)).fullid(10);
      unitnum = p(matching_files(3)).localid;
      [spikes3,epos3,evel3,names3,type3] = parsetrials(experiment,prefdirvec,unitnum,postwin);
  
      unit.allspikes = [spikes1;spikes2;spikes3];
      unit.allepos = [epos1;epos2;epos3];
      unit.allevel = [evel1;evel2;evel3];
      unit.allnames = [names1;names2;names3];
      unit.alltype = [type1 type2 type3];
  end
  
  unit.unitid = uniqueids(id);
  
  % find the preferred direction and set up for analysis of that direction

  for i = 1:8
    dex = find(unit.alltype == prefdirvec(i));
    unit.sp(i) = sum(sum(unit.allspikes(dex,:)));
  end
  
  unit.prefdir = prefdirvec(find(unit.sp == max(unit.sp)));
  
  % run the actual analysis, on the preferred direction, 
  % and save the image corresponding to the data
  dex = find(unit.alltype == unit.prefdir(1));
  unit.binaryspikes = unit.allspikes(dex,:); 
  unit.evel = unit.allevel(dex,:);

  % run on the de-saccaded eye velocity
  unit = tbtv_final(unit);
  print('-depsc2',[image_dir sprintf('unit%.3d.eps',uniqueids(id))]);
  close all
    
  % run the "forward" filter, from eye movement to spikes
  unit = tbtv_final_forward(unit);
  print('-depsc2',[image_dir sprintf('unit%.3d.forward.eps',uniqueids(id))]);
  close all

  % run on the eye position trace
  unit.epos = unit.allepos(dex,:);
  unit = tbtv_final_position(unit);

  % run on the aligned data
  unit = minimize_variance(unit);
  unit = tbtv_final_aligned(unit);
  
  % run for shorter data
  unit = tbtv_final_short(unit);
  
  % run to calculate the covariance
  
  % run on the flanks (off the preferred direction)
  % find the flanks above and below
  
  % flank above
  unit.abovedir = aboveflankvec(find(unit.sp == max(unit.sp)));
  dex = find(unit.alltype == unit.abovedir);
  unit.above.binaryspikes = unit.allspikes(dex,:); 
  unit.above.evel = unit.allevel(dex,:);
  unit = tbtv_final_above(unit);
  
  % flank below
  unit.belowdir = belowflankvec(find(unit.sp == max(unit.sp)));
  dex = find(unit.alltype == unit.belowdir);
  unit.below.binaryspikes = unit.allspikes(dex,:); 
  unit.below.evel = unit.allevel(dex,:);
  unit = tbtv_final_below(unit);
  
  % fit the preferred direction
  unit = fitprefdir(unit);
  
  % add the full ID
  unit.fullid = vertcat(p(matching_files).fullid);
  
  % and add the rank of the unit
  unit.rank = horzcat(p(matching_files).rank);
  
  % save the unit; change the variable's name to include the ID
  assignin('base',sprintf('unit%.3d',uniqueids(id)),unit)
  savefile = [output_dir sprintf('unit%.3d',uniqueids(id))];
  save(savefile,sprintf('unit%.3d',uniqueids(id)));
  clear(sprintf('unit%.3d',uniqueids(id)));
  
  disp(['Finished ID: ',num2str(uniqueids(id))]);
end
toc

%%%%%%%%%%%%
% Old Code %
%%%%%%%%%%%%

%  Here was code for winnowing the original set of ~300 cells down to ones
%  that actually responded appropriately for the study
%
%   if length(unit.prefdir) > 1
%     disp(['Unit ',num2str(uniqueids(id)),' has no single preferred direction']);
%     continue;
%   end
%   
%   % check for too few spikes
%   if max(unit.sp) < 25
%     disp(['Unit ',num2str(uniqueids(id)),' has fewer than 25 spikes']);
%     continue;
%   end
%     
%   % check for too few trials
%   if size(unit.evel,1) < 15
%     disp(['Unit ', num2str(uniqueids(id)),' has only ',num2str(size(unit.evel,1)), ' trials.'])
%     continue;
%   end
%   % check for too few spikes
%   if max(unit.sp) < size(unit.evel,1).*4
%     disp(['Unit ',num2str(uniqueids(id)),' has fewer than 4 spikes/trial, on average']);
%     continue;
%   end

  
% to run the analysis to accomodate the fact that trials were of
% different lengths, adapt the code below
%   if id <=  % unit 153 and below
%     postwin = 500;
%   elseif id >=  % unit 269 onward
%     postwin = 900;
%   else
%     postwin = 750; % between the two
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to create the pursuit database from the original unit database %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % load up the database of cells
% % originally from an excel file
% load unitdb;
% 
% % find the pursuit trials, and create a database of them
% % create a vector of unique ids, and of all ids.
% p = u(find(vertcat(u.experiment) == 'a'));
% tmp = find(vertcat(u.experiment) == 'a');
% for i = 1:length(tmp), p(i).unitdbind = tmp(i); end
%
% for i = 1:size(p,2)
%   cd(['/home/schoppik/exp/alldata/',p(i).fullid(1:8)]);
%   [names,types] = nex_info(strcat(p(i).fullid(1:end-7),'nex'));
%   os = min(find(types == 0)) - 1;
%   localid = strmatch(p(i).fullid(end-6:end),names) - os;
%   
%   if os < 1 | os > 12
%     disp('this should never have happened')
%   else
%     p(i).localid = localid;
%   end
% end

