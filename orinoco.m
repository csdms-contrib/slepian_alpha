function varargout=orinoco(res,buf)  
% XY=ORINOCO(res,buf)
% ORINOCO(...) % Only makes a plot
%
% Finds the coordinates of Orinoco, potentially buffered by some amount.
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
regn='orinoco';
c11=[260 15];
cmn=[345 -60];
xunt=[];

% This admittedly is a special case
XY=load(fullfile(getenv('IFILES'),'COASTS','orinoco'));

% Do it! Make it, load it, save it
XY=regselect(regn,XY(:,1),XY(:,2),xunt,res,buf);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

