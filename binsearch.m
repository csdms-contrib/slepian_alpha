function index=binsearch(x,sval)
% index=binsearch(x,sval)
%
% Binary search for unique values in sorted array.
%  
% INPUT:
%
% x      Sorted vector of numeric values
% sval   Numeric value to be searched in x
%  
% OUTPUT:
%
% index  index of sval with respect to x
%
% Last modified by fjsimons-at-alum.mit.edu, 18.08.2006
% (C) M. Khan mak2000sw@yahoo.com www.geocities.com/mak2000sw/ 
       
index=[];
from=1;
to=length(x);

while from<=to
    mid=round((from + to)/2);    
    dff=x(mid)-sval;
    if dff==0
        index=mid;
        return
    elseif dff<0 
        from=mid+1;
    else          
      to=mid-1;			
    end
end

