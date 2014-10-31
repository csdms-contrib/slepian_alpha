function defcell(name,value)
% function DEFCELL(name,value)
%
% A function similar to DEFVAL or DEFSTRUCT which will assign the default
% 'value' to a cell array 'name'. If the workspace or function calling 
% DEFCELL already has a cell variable 'name', it does not overwrite
% existing elements, but it will grow the existing cell to the size of
% the default, adding additional elements if necessary. If the named cell
% didn't exist at all, you can default it to anything you want in any
% number of dimensions. If it did exist, it will only work for up to
% two-dimensional arrays and complain otherwise.
%
% INPUT
%
% name    A string with a variable name
% value   The default cell array for that named variable
%
% OUTPUT:
%
%         None. The cell array named 'name' will "appear as if by magic
%            into your workspace or be available inside your function."
%
% See also ASSIGNIN, EVALIN, EXIST, ISEMPTY
%
% Last modified by gleggers-at-princeton.edu, 09/20/2013

% Error statement if a string is not provided as the first input (a very
% common mistake when typing quickly)
if ~ischar(name),
  error(sprintf(['The first argument of DEFCELL ',...
		 'has to be a string with a variable name']));
end

% Always do it is our default here
si=1;

% See if a cell named 'name' already exists in the calling function. If it
% does not, skip to the end and assign 'value' to a cell 'name' in the
% calling function. If it does, then there is a preexisting cell (an
% "original value" or "OV" as it will be abbreviated) that DEFCELL is
% trying to default over, and more careful consideration is necessary
% Remember, "exist" works on strings, "iscell" works on things!
if evalin('caller',[ 'exist(''' name ''',''var'')']) == 1 ...
  && evalin('caller',[ 'iscell(' name ')'])
  % Get size of the OV cell in the calling function
  ovS = evalin('caller',[ 'size(' name ')']);
  
  % Run through every element of that OV cell... (as a 2D  thing)
  for i = 1:ovS(1)
    for j = 1:ovS(2)
      % ...If an element (i,j) of that OV cell has a value, then 
      % replace the corresponding element (i,j) of this function's
      % 'value' cell with the element from the OV cell
      if evalin('caller',[ 'isempty(' name '{' num2str(i) ',' num2str(j) '})']) == 0
	value{i,j} = evalin('caller',[name '{' num2str(i) ',' num2str(j) '}']);
      end
    end
  end
else
  si=logical(0);
end
keyboard
% Assign this function's 'value' cell (modified or not) to variable 'name'
% in the calling function, replacing any existing variable of that name
if si
  assignin('caller',name,value);
end

