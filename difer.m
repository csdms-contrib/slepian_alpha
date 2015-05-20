function varargout=difer(V,tolex,sev,goods)
% DIFER(V,tolex,sev,goods)
% a=DIFER(V,tolex,sev,goods)
%
% Checks if the sum(abs(V(:))) exceeds 10^(-tolex)
%
% INPUT:
%
% V         Vector or matrix values (NOT symbolic ones!)
% tolex     Tolerance exponent (default: 10 for 1e-10)
% sev       0 produces WARNING upon failure (default)
%           1 produces ERROR upon failure 
% goods     A string with an uplifting message (if NaN, no message)
%
% OUTPUT (optional)
%
% a         1 if V is bigger than "zero", 0 if it isn't 
%           but remember if you request output there's no message 
%
% SEE ALSO: 
%
% ISEQUAL
%
% Last modified by fjsimons-at-alum.mit.edu, 06/27/2014

% Which one is the calling program?
[path,name,ext]=star69;

% Check in the calling routine that the incoming terms are not below
% machine precision even before they are being compared!
% disp('should make one called differr for relative comparisons')

defval('V',0)
defval('tolex',10)
defval('sev',0)
defval('goods',...
       sprintf('%s Check passed to tolerance E%i',...
	       upper(name),-tolex))
if ~isstruct(V) & ~strcmp(class(V),'sym')
  sabs=sum(abs(V(:)));
else
  error('No support for the struct/sym class')
end

% Make this explicitly equal to 1 if not symbolic variables will falsely test
if ~isnan(sabs)
  if sabs>10^(-tolex)
    mesg=sprintf('sum(abs(%s)) exceeds 0 by %8.3e',inputname(1),sabs);
    switch sev
     case 0
      if nargout==0
	warning(mesg)
      else
	varargout{1}=logical(1);
      end
     case 1
      if nargout==0
	error(mesg)
      else
	varargout{1}=logical(1);
      end
    end
  else
    if nargout~=0
      varargout{1}=logical(0);
    else
      if ~isnan(goods)
	disp(sprintf(goods,-tolex))
      end
    end
  end
else
  error('Comparison contains NaNs or is empty')
end
