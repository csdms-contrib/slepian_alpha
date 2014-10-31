function nrmorl=addmoff(nr,drk)
% nrmorl=ADDMOFF(nr,drk)
%
% Calculates the number of real spherical harmonic orders that belong to
% an expansion from degree l=0 to L, OR vice versa.
% For arrays where m=-l:l.
%
% INPUT:
%
% nr          The number of real spherical harmonic orders, OR
%             The degree of the expansion
% drk         'a' the input is the degree, calculate the degree [default]
%             'r' the input is the number, calculate the number
%
% OUTPUT:
%
% nrmorl      The degree of the expansion, OR
%             The number of real spherical harmonic orders
%
% EXAMPLE:
%
% In combination with GLMALPHA, use EL(addmoff(12)+1:addmoff(13))
% to find those coefficients of order 13.
%
% addmoff(m-1:L-1)+2*m+1 finds the location of the positive orders m in a
% vector listing the orders [0 0-11 0-11-22 0-11-22-33], as in, e.g. KERNELC
%
% See also: ADDMUP, ADDMOUT, ADDMON, ADDMIN
%
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011

defval('drk','a')

switch  drk
  case 'a'
    nrmorl=(nr+1).^2;
  case 'r'
    nrmorl=sqrt(nr)-1;
    if prod((round(nrmorl)==nrmorl)+0)==0
      warning('Invalid entry - noninteger result')
    end
 otherwise
  error('Specify a valid option')
end

