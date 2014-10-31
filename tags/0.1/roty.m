function Ry=roty(bita)
% Ry=ROTY(bita)
%
% Makes a rotation matrix for coordinates [x,y,z].
% Dahlen and Tromp, eq. (C.239).
%
% INPUT:
%
% bita     Rotation around y [radians]. If used on POINTS, 
%          the rotation is from positive x to positive z.
%
% OUTPUT:
% 
% Ry       The required rotation matrix
%
% EXAMPLE:
%
% newcords=[Ry*cords']'; % New coordinates of POINTS in old SYSTEM.
%
% EXAMPLE:
%
% x=0:10; y=repmat(0,1,length(x)); z=zeros(size(x)); clf ; plot3(x,y,z,'o-'); 
% bita=2*pi/6; xyzp=roty(bita)*[x ; y ; z]; hold on
% xp=xyzp(1,:)'; yp=xyzp(2,:)'; zp=xyzp(3,:)';
% plot3(xp,yp,zp,'rv-'); axis equal; axis([-10 10 -10 10 -10 10])
% xlabel('x'); ylabel('y'); zlabel('z'); grid on; view(130,30)
% title(sprintf('Rotation about y from blue to red over %i%s from +x to +z',...
%           round(bita*180/pi),str2mat(176)))
%
% SEE ALSO:
%
% ROTX, ROTZ, ROTS, ROTTP, ROTCOF, PLM2ROT
%
% Last modified by fjsimons-at-alum.mit.edu, 03/10/2010

Ry=[cos(bita) 0 -sin(bita);...
        0     1      0    ;...
    sin(bita) 0  cos(bita)];




