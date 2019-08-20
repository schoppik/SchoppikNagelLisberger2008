%
%
%
%
%          by: david schoppik
%        date: 08-15-2006
%     purpose: To go from the "unit???.mat" files to the "pair???.mat"
%              files
%
%   algorithm: find the number of unique days, calculate all possible
%              pairs, and then find out which correspond to significant
%              units


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make an experiment list %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

load pursdb
load unitdb
allids = vertcat(p.id);
expdb = [];

i = 2; % start here

cmp_exp = u(1).fullid(1:10);
next_exp = u(i).fullid(1:10);

expind = 1;
expdb(expind).numunits = 1;
expdb(expind).unitids = 1;

while i <= size(p,2)
  
  % is the next entry in the database from the same day?
  if strcmp(cmp_exp,next_exp)
    expdb(expind).id = next_exp;
    expdb(expind).numunits = expdb(expind).numunits + 1;
    expdb(expind).unitids(expdb(expind).numunits) = p(i).id;
    cmp_exp = next_exp;
    i = i + 1;
    if i < size(p,2)
      next_exp = p(i).fullid(1:10);
    end
    
  else % it was a new day
    expind = expind + 1;
    expdb(expind).id = next_exp;
    expdb(expind).numunits = 1; 
    expdb(expind).unitids(expdb(expind).numunits) = p(i).id;
    cmp_exp = next_exp;
    i = i + 1;
    if i <= size(p,2)
      next_exp = p(i).fullid(1:10);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a list of potential pairs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pp = []; 
for i = 1:size(expdb,2)
  if expdb(i).numunits > 2
    pp = horzcat(pp,expdb(i).unitids(generate_pairs(expdb(i).numunits)));
  elseif expdb(i).numunits == 2
    pp = horzcat(pp,expdb(i).unitids(generate_pairs(expdb(i).numunits))');
  end
end

% first, find the pairs that are repeats
% sometimes the same unit appears more than once during a single experiment
pp = pp(:,find(~ismember(1:size(pp,2),find(pp(1,:) == pp(2,:)))));

np = [];
% now find the pairs that are duplicates
% stupid nested loops
% use an "ismember" instead of an "and" to find pairs that are really the
% same with the order reversed
for i = 1:size(pp,2)
  for j = 1:size(pp,2)
   testpair = ismember(pp(:,i),pp(:,j));
   np(i,j) = and(testpair(1),testpair(2));
  end
end

% now eliminate the pairs that are dupes
for i = 1:size(np,2)
  keep(i) = min(find(np(i,:)));
end
keep = find(ismember(1:size(np,2),keep));

pairdb = pp(:,keep);

% now get rid of the pairs that correspond to units that didn't meet the
% initial criteria for analysis (i.e. ones that weren't significant)

load data
d = data.unitid;

pairdb = pairdb(:,find(ismember(pairdb(1,:),d)));
pairdb = pairdb(:,find(ismember(pairdb(2,:),d)));

goodpairdb = pairdb;

% one of the potential pairs is junk
% (it is a sorting/naming quirk)
pairdb = pairdb(:,[1:38,40:end]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the pair???.mat files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prefdirvec = 45:45:360;

% in case a pair is bad
offset = 0;

% now create the pair structures
for paircnt = 1:size(pairdb,2)
  pair = [];
  
  id1 = pairdb(1,paircnt);
  id2 = pairdb(2,paircnt);
  
  % load up the two units
  load(sprintf('unit%.3d',pairdb(1,paircnt)));
  load(sprintf('unit%.3d',pairdb(2,paircnt)));
  eval(sprintf('u1 = unit%.3d;',pairdb(1,paircnt)));
  eval(sprintf('u2 = unit%.3d;',pairdb(2,paircnt)));
  
  % find the trials in common
  
  % first, do the easy check: are all the trials the same?
  if strcmp(u1.allnames,u2.allnames)
    % sweet!
    pair = [];
    pair.spikes(:,:,1) = u1.allspikes;
    pair.spikes(:,:,2) = u2.allspikes;
    pair.epos = u1.allepos;
    pair.evel = u1.allevel;
    pair.type = u1.alltype;
    pair.spikesbydir(:,1) = u1.sp;
    pair.spikesbydir(:,2) = u2.sp;
    pair.prefdir(:,1) = u1.prefdir;
    pair.prefdir(:,2) = u2.prefdir;
    pair.id = pairdb(:,paircnt);
         
    mf1 = find(allids == pair.id(1));
    mf2 = find(allids == pair.id(2));
    
    mf1 = mf1(1);
    mf2 = mf2(1);
    
    e1 = str2num(p(mf1).fullid(end-1));
    e2 = str2num(p(mf2).fullid(end-1));
    pair.electrodedist = abs(e1-e2);
    
  else
    % now we have to figure out the indices of the trials that are in
    % common
    
    % first, figure out which one has fewer trials
    [jnk smaller] = min([size(u1.allnames,1),size(u2.allnames,1)]);
    dex = [];
    switch smaller
      case 1 % u1 has fewer trials
        % make a cell array of the trial names for the longer one
        names = cellstr(u2.allnames);
        for namecnt = 1:size(u1.allnames,1)
          match = find(strcmp(u1.allnames(namecnt,:),names)); 
          if ~isempty(match)
            dex(namecnt) = match;
          end
        end
       u1dex = find(dex);
       u2dex = dex(find(dex));
           
      case 2
        % make a cell array of the trial names for the longer one
        names = cellstr(u1.allnames);
        for namecnt = 1:size(u2.allnames,1)
          match = find(strcmp(u2.allnames(namecnt,:),names));
          if ~isempty(match)
            dex(namecnt) = match;
          end
        end
       u2dex = find(dex);
       u1dex = dex(find(dex));
    end

    pair = [];

    pair.spikes(:,:,1) = u1.allspikes(u1dex,:);
    pair.spikes(:,:,2) = u2.allspikes(u2dex,:);
    pair.epos = u1.allepos(u1dex,:);
    pair.evel = u1.allevel(u1dex,:);
    pair.type = u1.alltype(u1dex);
    
    pair.prefdir(:,1) = u1.prefdir;
    pair.prefdir(:,2) = u2.prefdir;
    pair.id = pairdb(:,paircnt);
    
    mf1 = find(allids == pair.id(1));
    mf2 = find(allids == pair.id(2));
    
    mf1 = mf1(1);
    mf2 = mf2(1);
    
    e1 = str2num(p(mf1).fullid(end-1));
    e2 = str2num(p(mf2).fullid(end-1));
    pair.electrodedist = abs(e1-e2);
    
  end

  
  % if we're done, run the nbnv analysis
  pair = nbnv_final_tbtv(pair);
  
  
  % save the pair
  assignin('base',sprintf('pair%.3d',paircnt),pair);
  savefile = sprintf('/home/schoppik/exp/long/pair%.3d',paircnt);
  save(savefile,sprintf('pair%.3d',paircnt));
  
  % clear out the old stuff
  clear(sprintf('pair%.3d',paircnt));
  clear(sprintf('unit%.3d',pairdb(1,paircnt)));
  clear(sprintf('unit%.3d',pairdb(2,paircnt)));
  clear('u1');
  clear('u2');
  
  disp(paircnt)
  
end % the big loop
  
%   if size(pair.spikes,1) < 120
%     disp(sprintf('too few trials in pair %.3d',paircnt))
%     goodpairdb = goodpairdb(:,[1:paircnt-1 paircnt+1:end]);
%     continue;
%   end

