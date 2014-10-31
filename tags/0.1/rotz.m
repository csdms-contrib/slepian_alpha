function Rz=rotz(gama)
% Rz=ROTZ(gama)
%
% Makes a rotation matrix for coordinates [x,y,z].
% Dahlen and Tromp, eq. (C.240).
%
% INPUT:
%
% gama     Rotation around z [radians]. If used on POINTS,
%          the rotation is from positive y to positive x.
%
% OUTPUT:
% 
% Rz       The required rotation matrix
%
% EXAMPLE:
%
% newcords=[Rz*cords']'; % New coordinates of POINTS in old SYSTEM. 
%
% EXAMPLE:
%
% x=0:10; y=0:10; z=zeros(size(x)); clf; plot3(x,y,z,'o-'); 
% bita=pi/6; xyzp=rotz(bita)*[x ; y ; z]; hold on
% xp=xyzp(1,:)'; yp=xyzp(2,:)'; zp=xyzp(3,:)'; 
% plot3(xp,yp,zp,'rv-'); axis equal; axis([-10 10 -10 10 -10 10])
% xlabel('x'); ylabel('y'); zlabel('z'); grid on; view(130,30)
% title(sprintf('Rotation about z from blue to red over %i%s from +y to +x',...
%           round(bita*180/pi),str2mat(176)))
%
% EXAMPLE: 
%
% Illustrate the trick to make the polar rotation over +/-90 only,
% which, for spherical harmonic rotation, makes all the difference....
% 
% a=randn; b=randn; g=randn; 
% difer(rotz(g)*roty(b)*rotz(a)-...
%       rotz(g+pi/2)*roty(pi/2)*rotz(0)*rotz(b)*roty(-pi/2)*rotz(a-pi/2))
%
% SEE ALSO:
%
% ROTX, ROTY, ROTCOF, DLMB, BLANCO, PLM2ROT
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/10/2010

Rz=[ cos(gama) sin(gama)  0;...
    -sin(gama) cos(gama)  0;...
        0          0      1];


