function varargout=honshu(res,buf)
% XY=HONSHU(res,buf)
% HONSHU(...)  % Only makes a plot
%
% Finds the coordinates of Honshu and Kyushu, potentially buffered by
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
% Last modified by fjsimons-at-alum.mit.edu, 09/24/2014

defval('res',0)
defval('buf',0)

% Parameters that make this the region in question
regn='honshu';
c11=[129 42];
cmn=[146 30];
xunt=[3:56];
ofs=0;

% Do it! Make it, load it, save it
XY=regselect(regn,c11,cmn,xunt,res,buf,ofs);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

