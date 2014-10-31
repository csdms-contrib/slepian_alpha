function varargout=africa(res,buf) 
% XY=AFRICA(res,buf)
% AFRICA(...) % Only makes a plot
%
% Finds the coordinates of Africa, potentially buffered by some amount.
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
% Last modified by charig-at-princeton.edu, 11/23/2011
% Last modified by fjsimons-at-alum.mit.edu, 11/23/2011

defval('res',0)
defval('buf',0)

% Parameters that make this the region in question
regn='africa';
c11=[0 37.5 ; 342.5 36.5];
cmn=[52 -35 ; 360 3.5];
xunt=[261:307 110:260 339:399];
ofs=[360 0];

% Do it! Make it, load it, save it
XY=regselect(regn,c11,cmn,xunt,res,buf,ofs);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);
