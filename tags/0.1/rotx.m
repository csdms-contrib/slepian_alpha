function Rx=rotx(alfa)
% Rx=ROTX(alfa)
%
% Makes a rotation matrix for coordinates [x,y,z].
% NOT one of Dahlen and Tromp's (C.238-240).
%
% INPUT:
%
% alfa     Rotation around x [radians]. If used on POINTS,
%          the rotation is from positive z to positive y.
%
% OUTPUT:
% 
% Rx       The required rotation matrix
%
% EXAMPLE:
%
% newcords=[Rx*cords']'; % New coordinates of POINTS in old SYSTEM.
%
% EXAMPLE:
%
% x=repmat(0,1,length(x)); y=0:10; z=0:10; clf; plot3(x,y,z,'o-'); 
% bita=pi/6; xyzp=rotx(bita)*[x ; y ; z]; hold on
% xp=xyzp(1,:)'; yp=xyzp(2,:)'; zp=xyzp(3,:)'; 
% plot3(xp,yp,zp,'rv-'); axis equal; axis([-10 10 -10 10 -10 10])
% xlabel('x'); ylabel('y'); zlabel('z'); grid on; view(130,30)
% title(sprintf('Rotation about x from blue to red over %i%s from +z to +y',...
%           round(bita*180/pi),str2mat(176)))
%
% SEE ALSO:
%
% ROTX, ROTZ, ROTS, ROTTP, ROTCOF, PLM2ROT
%
% Last modified by fjsimons-at-alum.mit.edu, 03/10/2010

Rx=[1     0          0    ;...
    0  cos(alfa) sin(alfa);...
    0 -sin(alfa) cos(alfa)];

