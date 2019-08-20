% generate a covariance matrix
%
%
%
%         by: david schoppik
%       date: 4/23/07
%    purpose: to generate appropriate covariance matrices with known
%             off-diagonal mean and variance
%
%      usage: [C m v] = smartCov(n,eij,vij,tol);
%
%             C is an n-by-n matrix with mean m and variance v where:
%             eij is the off-diagonal mean
%             vij is the off-diagonal variance
%             tol is the allowable error - if not present, doesn't check
%
%
%  based on formulae 6,7,19, in Hirschberger et. al. "Randomly Generating
%  Portfolio-Selection Covariance Matrices with Specified Distributional
%  Characteristics," this function calculates a matrix F such that F*F' has
%  the desired off-diagonal mean and variance
%
%

function [C m v] = smartCov(n,eij,vij,tol);

% check the inputs
if nargin < 3
  help smartCov
  return;
end

if nargin == 3;
  tol = 0;
end

if eij < 0 | vij < 0 |  n < 0
  error('All quantities must be > 0')
end

% variables
maxloops = 10000;
numloops = 0;
eii = 1;
verbose = 0;
returnCorr = 1; % return the correlation matrix, not the covariance

% calculate the desired features of the distribution that will make up
% matrix F

% formula 19
m = round( (eii.^2 - eij.^2) ./ vij );

% formula 6
ehat = sqrt(eij/m);

% formula 7
vhat = -eij.^2 + sqrt(eij.^4 + vij/m);

% generate and transform a set of draws into a covariance matrix

% use a standard normal
qij = randn(n,m);
fij = ehat + (sqrt(vhat)*qij);

% use a log-normal distribution
%fij = lognrnd(ehat,vhat,m,n);

C = fij*fij';

% if a tolerance has been specified
if tol ~= 0
  jnk = tril(C,-1);
  jnk = jnk(find(jnk));

  while abs(mean(jnk)-eij) > tol | abs(var(jnk)-vij) > tol

    % draw again
    % use a standard normal
    qij = randn(n,m);
    fij = ehat + (sqrt(vhat)*qij);
    
    % use a log-normal distribution
    % fij = lognrnd(ehat,vhat,m,n);

    C = fij*fij';

    jnk = tril(C,-1);
    jnk = jnk(find(jnk));

    numloops = numloops + 1;
    if numloops > maxloops
      disp('Could not match requested criteria')
      C = nan;
      return;
    else
    end
 
  end
end


if returnCorr == 1;
  C = diag(1./sqrt(diag(C))) * C * diag(1./sqrt(diag(C)));
end

if verbose == 1;
  jnk = tril(C,-1);
  disp(sprintf('Off-diagonal mean = %0.4g, variance = %0.4g',...
    mean(jnk(find(jnk))),var(jnk(find(jnk)))))
end

if nargout == 3
   jnk = tril(C,-1);
   m = mean(jnk(find(jnk)));
   v = var(jnk(find(jnk)));
end
