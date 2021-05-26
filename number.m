function varargout=number(strmat)
% numed=NUMBER(strmat)
%
% Generates/displays row numbers to index a (string) matrix or cell array
%
% INPUT:
%
% strmat         a certain input matrix or cell array
%
% OUTPUT:
%
% numed          a string matrix with the row numbers
%                if no output requested, displays to screen
%
% EXAMPLES:
%
% number([1 2 3; 3 4 5])
% number([ 'menu' ; 'opti'   ])
% number([{'menu'} {'option'}])
% number({{'menu'} {'option'}})
%
% fjsimons-at-alum.mit.edu, 05/26/2021

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

% If no output, just displays
if nargout==1
  varargout{1}=numed;  
else
  disp(numed)
end

