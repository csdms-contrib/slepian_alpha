function na=axisc()
%AXISC Switch back to last axis coordinates.
%
%   H = AXISC sets as current and returns the handle of the axis that was
%   current before the call to FIGC.
% 
%   See also FIGC

% Written by: Oded Aharonson, 10/28/1999
global oldaxisglobal

if length(oldaxisglobal)>0
  axes(oldaxisglobal);
  if (nargout > 0) 
    na=oldaxisglobal;
  end
else
  if (nargout > 0) 
    na=axes;
  else
    axes;
  end
end
