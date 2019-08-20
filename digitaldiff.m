function dxdt = digitaldiff(x,dt,lag)
%
%
%       by: John G. O'Leary & David Schoppik
%     date: March 10, 2005
%  purpose: A zero-delay digital differentiator
%
%
%  dxdt = (x(t+lag) - x(t-lag)) ./ ((2*lag)*dt)
%
% defaults: lag = 1 sample, dt = .001 sec
switch nargin 
  case {2}
  lag = 1;

  case {1}
  lag = 1;
  dt = .001;
  
  case {0}
    help digitaldiff
    return
end

if lag < 1
  disp('lag has to be an integer > 1')
  return
end


if length(x) < (2*lag) 
  disp('input x is too short -- must be longer than 2*lag')
  return
end

% force lag into an integer
lag = round(lag);

dxdt = [zeros(1,lag) x(1+2*lag:end) - x(1:end-2*lag) zeros(1,lag)] ./ ((2*lag)*dt);

