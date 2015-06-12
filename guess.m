function guesses=guess(n)
% GUESS(n)
%
% Returns integer guesses between -10 and +10 including 0
%
% INPUT:
%
% n            How many guess you'd like to get [default: 10]
%
% OUTPUT:
%
% guesses      A bunch of n random integers
%
% SEE ALSO: 
%
% RANDIJ
%
% Last modified by fjsimons-at-alum.mit.edu, 2/6/2009

defval('n',10)

guesses=ceil((-1).^ceil(randn(n,1)).*(10*rand(n,1)))';
