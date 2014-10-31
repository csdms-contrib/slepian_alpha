function varargout=plotnorm(x,varargin)
% h=plotnorm(x,varargin)
%
% Plots data on length interval [0 1]

if ~nargin
  h=plot(linspace(0,1,length(x)),x);
else
  h=plot(linspace(0,1,length(x)),x,varargin{:});
end

if nargout==1
  varargout{1}=h;
end
