function varargout=baffin(res,buf)  
% XY=baffin(res,buf)
% baffin(...) % Only makes a plot
%
% Finds the coordinates of Baffin Island, potentially buffered by some amount.
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% buf      Distance in degrees that the region outline will be enlarged
%          by BUFFERM, not necessarily integer, possibly negative
%          [default: 0]
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
%
% Last modified by fjsimons-at-alum.mit.edu, 09/23/2014

defval('res',0)
defval('buf',0)

% Parameters that make this the region in question
regn='baffin';
c11=[270 77];
cmn=[300 58];
xunt=[14:107];

% Do it! Make it, load it, save it
XY=regselect(regn,c11,cmn,xunt,res,buf);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);
