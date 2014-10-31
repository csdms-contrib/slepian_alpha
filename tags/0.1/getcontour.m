function [bigc,val]=getcontour(c)
% [bigc,val]=GETCONTOUR(c)
%
% From a contour matrix, gets the contours
%
% INPUT:
%
% c       A contour matrix out of CONTOUR or CONTOURC
% 
% OUTPUT:
%
% bigc    A cell array with the coordinates of these contours
% val     A vector with the contour values
%
% EXAMPLE:
%
% for index=1:length(val); plot(bigc{index}(1,:),bigc{index}(2,:))
% hold on; end
%
% Last modified by fjsimons-at-alum.mit.edu, 04/20/2010

nxy=size(c,2);

b=2;
index=0; e=0;
while e<nxy
  index=index+1;
  val(index)=c(1,b-1);
  e=b+c(2,b-1)-1;
  bigc{index}=c(1:2,b:e);
  b=e+2;
end

