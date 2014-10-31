function na=figc()
%FIGC Switch to figure normalized coordinates
%
%   H = FIGC sets as current and returns the handle of an invisible axis
%   with normalized figure coordinates.
% 
%   See also AXISC

% Written by: Oded Aharonson, 10/28/1999

global oldaxisglobal
if (length(get(gcf,'children'))>0) 
  oldaxisglobal=gca;
else
  oldaxisglobal=[];
end

nna=axes('position',[0 0 1 1], ...
    'color','none','xtick',[],'ytick',[],'visible','off');

if (nargout > 0) 
  na=nna;
end