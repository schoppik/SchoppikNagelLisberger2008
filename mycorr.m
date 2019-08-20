function r = mycorr(A,B);
% not expecting things to be flipped
% i.e. time is across rows
% will only return the correlation along the diagonal
% which is identical to the output of corr;

if nargin == 2;
% get rid of the means
A = A-repmat(mean(A),size(A,1),1);
B = B-repmat(mean(B),size(B,1),1);

% compute the covariance
num = diag((A'*B));
denom = sqrt(diag(A'*A).*diag(B'*B));
warning off
r = num./denom;
warning on

else
  B = A(:,2);
  A = A(:,1);
  
  A = A-repmat(mean(A),size(A,1),1);
  B = B-repmat(mean(B),size(B,1),1);
  
  % compute the covariance
  num = diag((A'*B));
  denom = sqrt(diag(A'*A).*diag(B'*B));
  warning off
  r = num./denom;
  warning on
end