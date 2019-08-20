%
%
%
%            by: david schoppi
%          date: 8-8-06
%       purpose: to minimize the variance between a matrix of eye movements
%
%        usage: unit = minimize_variance(unit) where unit is the standard
%               unit structure

function unit = minimize_variance(unit);

range = 100;

% first, shuffle the cells, and choose a seed
dex = randperm(size(unit.evel,1));
seed = unit.evel(dex(1),:);
unit.shift(dex(1)) = 0;


for j = 2:length(dex)
  sse = [];


  for i = -range:range
    % reassign the value so the circshift is the right amount
    eval = unit.evel(dex(j),:);
    eval = circshift(eval,[0 i]);

    % take the sum of the squared difference, calculated both ways
    % (cell-eval) and (eval-cell)
    sse(i+range+1) = sum((eval(range+1:end-range) - seed(range+1:end-range)).^2) + ...
      sum((seed(range+1:end-range) - eval(range+1:end-range)).^2);
  end

  
  [jnk shift] = min(sse);
  shift = shift - (range + 1);
  
  if shift < range & shift > -range % in case it didn't converge
    unit.shift(dex(j)) = shift;
  else
    unit.shift(dex(j)) = 0;
    disp('did not converge')
  end
  
  % now recalculate the seed as the mean of the cells that have been run thus
  % far

  eval = unit.evel(dex(j),:);
  seed = (((j-1).*seed) + circshift(eval,[0 unit.shift(dex(j))]))./j;

end
