function f=norm3d(x)
% f=NORM3D(x)
%
% Calculates the 2-norm of a three-dimensional array in x
%
% Last modified by Ignace Loris (igloris@vub.ac.be) on 22.06.2009
% Last modified by fjsimons-at-alum.mit.edu on 06/30/2009

f=sqrt(sum(sum(sum(x.^2))));
