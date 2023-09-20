function varargout=ukraine(res,buf)
% XY=UKRAINE(res,buf)
% UKRAINE(...) % Only makes a plot
%
% Finds the coordinates of Ukraine, potentially buffered by some amount.
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
% Last modified by fjsimons-at-alum.mit.edu, 09/20/2023

defval('res',0)
defval('buf',0)

% Parameters that make this the region in question
regn=mfilename;
xunt=[];

% This admittedly is a special preloaded case
XY=load(fullfile(getenv('IFILES'),'COASTS',regn));

% Modify and resave it
XY=regselect(regn,XY(:,1),XY(:,2),xunt,res,buf);

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

