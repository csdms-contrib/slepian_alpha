function ig=trapeze(z,f)
% ig=TRAPEZE(z,f)
%
% Quadrature of a function using the trapezoidal rule.
%
% INPUT:
%
% z     Argument (vector)
% f     Function values (vector or matrix)
%
% OUTPUT:
%
% ig    Integral of f bewteen z(1) and z(end) 
%
% SEE ALSO: SIMPSON
%
% Last modified by fjsimons-at-alum.mit.edu, 02/22/2007

if size(z,1)==1; z=z(:); end
if size(f,1)==1; f=f(:); end

[m,n]=size(f);

ig=nan(n,1);

for index=1:n
  ig(index)=abs(diff(z))'*(f(1:m-1,index)+f(2:m,index))/2;
end

