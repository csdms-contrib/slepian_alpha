function [l,m,mu,check,tol]=pxyerh(l,m,mu,check,tol)
% [l,m,mu,check,tol]=PXYERH(l,m,mu,check,tol)
%
% Error handling for PLM, XLM, YLM.
%
% Not a stand-alone program. See PLM, XLM, and YLM.
%
% Last modified by fjsimons-at-alum.mit.edu, 03/16/2016

% Input error catching
if (prod(size(m))~=1 | isempty(m)) & prod(size(l))~=1
  error('Degree and order cannot both be vectors')
end
if any(mu>1 | mu<-1)
  error('Argument exceeds allowable range')
end
if any(l)<0
  error('Order must be positive')
end
if any(abs(m)>abs(l))
  error('Order exceeds allowable range')
end
if isempty(m)
  m=0:l;  
end
% Must be a row vector for the zero-degree case
if min(size(mu))==1
  mu=mu(:)';
end
