function x=indeks(M,i)
% x=INDEKS(M,i)
%
% Extracts linearly indexed position(s) i out of a matrix M.
% Works for logical and numeric indices.
%
% INPUT:
%
% M         The input matrix 
% i         The requested set of running linear indices [default: 1]
%
% EXAMPLES:
% 
% indeks([1 2],2) 
% indeks([1 2],':,2')
% indeks([1 2],'end')
%
% See also RINDEKS, KINDEKS, TINDEKS, DINDEKS, SQUEEZE
%
% Last modified by fjsimons-at-alum.mit.edu, 01/30/2019

defval('i',1)

if ~isstr(i)
  x=M(i(:));
else
  eval([ 'x=M(' i ');'])
end
