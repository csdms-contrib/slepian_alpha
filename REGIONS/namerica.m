function varargout=namerica(res,buf)  
% XY=NAMERICA(res,buf)
% NAMERICA(...) % Only makes a plot
%
% Finds the coordinates of Namerica, potentially buffered by some amount.
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
regn='namerica';
c11=[191   74.5];
cmn=[304.5 12  ];
xunt=[385:717 721:1004];

% Do it! Make it, load it, save it
XY=regselect(regn,c11,cmn,xunt,res,buf);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);
