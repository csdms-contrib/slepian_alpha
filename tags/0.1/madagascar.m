function varargout=madagascar(res,scal)
% XY=MADAGASCAR(res,scal)
% MADAGASCAR(...)  % Only makes a plot
%
% Finds the coordinates of Madagascar, potentially buffered by
% some amount. 
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
% XY       Closed-curved coordinates of the island
%
% Last modified by fjsimons-at-alum.mit.edu, September 28, 2014

defval('res',0)
defval('buf',0)

% Parameters that make this the region in question
regn='madagascar';
c11=[41 -11];
cmn=[52 -27];
xunt=[1:35];
ofs=0;

% Do it! Make it, load it, save it
XY=regselect(regn,c11,cmn,xunt,res,buf,ofs);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

