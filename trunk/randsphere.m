function varargout=randsphere(N,D,rad)
% [x,y,z]=RANDSPHERE(N,D,rad)
% [lon,lat]=RANDSPHERE(N,D,rad)
%
% Generates random points on a D-dimensional unit sphere,
% either in Cartesian or spherical coordinates.
%
% INPUT:
%
% N         Desired number of random points [default: 100]
% D         Dimension [default: 3, for the sphere]
% rad       Radius of the sphere [default: 1]
%
% OUTPUT:
%
% x,y,z     Cartesian points on the unit sphere, OR
% lon,lat   Earth coordinates (in degrees)
%
% Last modified by fjsimons-at-alum.mit.edu, 04/20/2009

defval('N',100)
defval('D',3)
defval('rad',1)

if nargout==2 & D==3
  % Pick random points on a sphere
  % Not just phi uniform on 0->2pi and theta uniform on 0->pi
  % If that was the case, the probability of falling in patch
  % tended by dtheta by dphi is sin(theta)dtheta*dphi
  % (infinitesimal solid angle, area element...).
  % So the smaller theta, the larger the probability, which leads
  % to a clustering at the poles.
  % But if pick theta=acos(1-2u) where u is uniform on 0->1
  % then dtheta=-1/(sqrt(1-u^2))du, 
  % sin(theta)=sin(acos(1-2u))=sqrt(1-u^2) so
  % sin(theta)*dthat=-du
  % hence probability of falling in this patch is -du*dphi/4pi
  % which is clearly uniform.
  U=rand(N,1);
  V=rand(N,1);
  ph=2*pi*U;
  th=acos(2*V-1);
  varargout{1}=ph*180/pi;
  varargout{2}=(pi/2-th)*180/pi;
elseif nargout==3 | D~=3
  % Number of dimensions; generalized for higher dimensions
  % Marsaglia, G. "Choosing a Point from the Surface of a Sphere"
  % Ann. Math. Stat. 43, 645-646, 1972.
  % Muller, M. E. "A Note on a Method for Generating Points Uniformly 
  % on N-Dimensional Spheres" 
  % Comm. Assoc. Comput. Mach. 2, 19-20, 1959. 
  R=randn(N,D);         
  R=diag(1./sqrt(sum(R'.^2)))*R;
  for index=1:D
    varargout{index}=R(:,index)*rad;
  end
end



