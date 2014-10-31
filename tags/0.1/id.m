function varargout=id(str,varargin)
% ID(str)
% ID(str,color)
% handle=ID(...)
%
% Function to place name, date, and calling function name on plot bottom.
% If 'str'=[] then default copyright marker.
%
% EXAMPLE:
%
% id;
%
% or
% id('bad plot!');
%
% Last modified by fjsimons-at-alum.mit.edu, 10/17/2011

defval('str',[])
if isempty(str)
  str=copyright;
  [path,prog,ext]=star69;
end

figc
if strcmp(orient,'landscape')
  a=text(0.9412,0.02,str,'vert','bottom','horiz','right','fontsize',9);
elseif strcmp(orient,'tall')
  a=text(0.9195,0.0455,str,'vert','bottom','horiz','right','fontsize',9);  
elseif strcmp(orient,'portrait')
  a=text(0.95,0.025,str,'vert','bottom','horiz','right','fontsize',9);  
end
if nargin>1
  set(a,'Color',varargin{1})
end
axisc

% Output
varns={a};
varargout=varns(1:nargout);
