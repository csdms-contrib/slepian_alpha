function varargout=number(strmat)
% number(strmat)
% numed=number(strmat)
%
% "Numbers" or indexes a string matrix or cell by the lines
%
% Uses a recursive algorithm, kind of
%
% fjsimons-at-alum.mit.edu, 12/17/02

if ~iscell(strmat)
  if ~ischar(strmat)
    strmat=int2str(strmat);
  end
  [m,n]=size(strmat);
  nums=num2str((1:m)');
  numed=[nums repmat(str2mat(32),m,1) strmat];
elseif iscell(strmat)
  numed=number(strcell(strmat));
end

if nargout==1
  varargout{1}=numed;  
else
  disp(numed)
end

