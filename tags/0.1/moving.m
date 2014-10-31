function varargout=moving(y,wlen)
% [ma,mo]=MOVING(y,wlen)
% 
% Moving average routine of signal 'y' with window length 'wlen'
%
% INPUT:
%
% y       The signal (which gets vectorized)
% wlen    The moving window length (in samples)
%
% OUTPUT:
%
% ma      The signal after moving averaging (not original length!)
% mo      The missing piece to restore he length, simply the repeated last entry
%
% Written by fjsimons-at-alum.mit.edu, 03/06/2014

y=y(:);
m=length(y);

if m < wlen
  error('Window longer than the signal itself');
else
  % Don't attempt to fix for NaN's
  ma=cumsum([sum(y(1:wlen)); y(wlen+1:m)-y(1:m-wlen)])./wlen;
  if nargout>1
    mo=repmat(ma(end),wlen-1,1);
  else
    mo=[];
  end
end


% Optional arguments
varns={ma,mo};
varargout=varns(1:nargout);
