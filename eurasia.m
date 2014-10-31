function varargout=eurasia(res,buf)  
% XY=EURASIA(res,buf)
% EURASIA(...) % Only makes a plot
%
% Finds the coordinates of Eurasia, potentially buffered by some amount.
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
regn='eurasia';
c11=[0 77.5 ; 350 50 ; 180 69];
cmn=[180 8; 360 36 ; 190 64];
xunt=[420:827 914:1023 1556:1605 1030:1548 1607:1639];
ofs=[360 0 360];

% Do it! Make it, load it, save it
XY=regselect(regn,c11,cmn,xunt,res,buf,ofs);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

