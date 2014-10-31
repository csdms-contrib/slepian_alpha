function I=simpson(z,f)
% I=SIMPSON(z,f)
%
% Numerical integration
%
% INPUT:
%
% z     Arguments of the function f, can be unequally spaced, or can be
%       in pairs of intervals with the spacing varying by pair
% f     The function values f(z), can be a matrix with functions down the columns
%
% OUTPUT:
%
% I     The value of the integral between z(1) end z(end)
%       computed by Simpson's (if z is equally spaced) or trapezoidal
%       rule (if z is unequally spaced)
%
% EXAMPLES:
%
% z=[0 1 2 4 6 7 8 13 18]'; f=[sin(z) cos(z)];
% I=simpson(z,f)
% z=[0 1 2 4 6 7 8.02 13 18]'; f=[sin(z) cos(z)];
% I=simpson(z,f)
% z=linspace(0,18,32)'; f=[sin(z) cos(z)];
% I=simpson(z,f)
% z=[1 5 6 9 11 14]'; f=[sin(z) cos(z)];
% I=simpson(z,f)
%
% Last modified by fjsimons-at-alum.mit.edu, 10/26/2011

z=z(:);
if size(f,1)==1; f=f(:); end

[m,n]=size(f);

% Odd is the sequence of choice here
% Note that 'z' defines PAIRS of layers with equal thickness (as surfc
% models usually have). 
% It checks for that but allows for two variations
% of the equal thickness that then must not exceed
% the next bigger interval/100.
zd=diff(diff(z)); 
zd=abs(zd(1:2:end));  
zd=zd(~~zd);

if ~isempty(zd)
  uzd=unique(diff(z)); uzd=uzd(3);
  if sum(zd)
    if any(zd>=repmat(uzd/100,length(zd),1))
      disp('Reduced order of integration method to trapezoidal rule')
      I=trapeze(z,f);
      return
    end
  end
end 

% Note Numerical Recipes 4.1.4 vs 4.1.13... alternating 1 4 2 4 2 4 1
% on successive non-overlapping pairs of intervals...
for index=1:n
  I(index)=[z(3:2:m)-z(1:2:m-2)]'*...
      [f(1:2:m-2,index)+4*f(2:2:m-1,index)+f(3:2:m,index)]/6;
end
