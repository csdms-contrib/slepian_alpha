function w=rotenberg(e,s)
% Converts Rotenberg's peculiar notation to a numerical value
%
% INPUT:
%
% e    The vector with exponents on the prime numbers
% s    0 positive sign
%      1 negative sign
%
% OUTPUT:
%
% w    The numerical value of these symbols
%
% Last modified by fjsimons-at-alum.mit.edu, 12/28/2006

defval('s',0)

P=primes(100);

w=(-1)^s*sqrt(prod(P(1:length(e)).^e));



