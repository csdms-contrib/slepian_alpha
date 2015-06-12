function varargout=rottp(theta,phi,alfa,bita,gama)
% [thetap,phip,rotmats]=ROTTP(theta,phi,alfa,bita,gama)
%
% Rotates POINTS on the unit sphere surface using Euler angles in an
% ACTIVE rotation convention (the inverse of DT pp 920-924). The results
% are coordinates for the new points, in the original, unrotated
% coordinate system. To invert, flip signs and orders of the Euler angles.
%
% INPUT:
%
% theta            Colatitudes of the set of points [radians]
% phi              Longitudes of the set of points [radians]
%                  If theta and phi are not of the same size, theta is
%                  interpreted as a matrix with constant rows, and phi
%                  as a matrix with constant columns.
% alfa, bita, gama Euler angles [radians], POINTS are rotated over 
%                     alf (0<2pi) around z increasing from y to x, then 
%                     bta (0<=pi) around old y increasing from x to z, then
%                     gam (0<2pi) around old z increasing from y to x.
%                  Equivalently but alternatively, POINTS are rotated over
%                     gama (0<2pi) around z increasing from y to x, then
%                     bita (0<=pi) around new y increasing from x to z, then
%                     alfa (0<2pi) around new z increasing from y to x.
%                  Or equivalently, rotates POINTS over
%                     gama-pi/2 around new z increasing from y to x, then
%                          pi/2 around new y increasing from x to z, then
%                          0    around new z increasing from y to x.
%                  followed by a second rotation of the POINTS,
%                     beta      around z increasing from y to x, then
%                         -pi/2 around new y increasing from x to z, then
%                     alfa+pi/2 around new z increasing from y to x, then
%                  Note that in order to do perform this ACTIVE rotation
%                  as described above, we stick to the PASSIVE rotation
%                  framework described by Dahlen and Tromp, App. C.8., in
%                  other words the equations for PASSIVE are the ones for
%                  ACTIVE but backwards, flipping their signs and
%                  orderings. Hence the perversity of writing about
%                  rotating points and finding the expansion in the
%                  unrotated system but talking about "new" axes in the
%                  description... it's AS IF the original axes were
%                  rotated like the "points", which we use strictly to
%                  clarify the operations...
%                  Find the appropriate rotation angles for simple
%                  geographic location, that is, the former North Pole moves:
%                  alp=0;      % Around original z axis, "clockwise"
%                  bta=lat-90; % To desired latitude, around old y axis
%                  gam=180-lon;   % To desired longitude, around old z axis 
%                  which should be consistent with the arguments in PTOSLEP
%
% OUTPUT:
%
% thetap           New colatitudes of the points [radians] in the old system
% phip             New longitudes of the points [radians] in the old system 
% rotmats          The actual rotation matrix used in the transformation
%
% SEE ALSO: PLM2ROT
%
% SEE: Dahlen and Tromp (DT), Appendix C.8.
%
% EXAMPLE:
%
% rottp('demo1')
% rottp('demo2')
% rottp('demo3')
% 
% EXAMPLE:
%
% clf; [a,b,lola]=plotcont([90 -10.7],[160 -40]); 
% lola=lola(~isnan(lola(:,1)) & ~isnan(lola(:,2)),:)*pi/180;
% angs=linspace(0,2,10);
% for ix=angs
% [colap,lop]=rottp(pi/2-lola(:,2),lola(:,1),0,ix,ix);
% hold on; plot(lop*180/pi,(pi/2-colap)*180/pi,'b');
% end; hold off
%
% theta=linspace(pi/3,pi/2,30); phi=linspace(0,2*pi/3,31);
% [thetap,phip]=rottp(theta,phi,1/2,1/3,1/4);
% [phi,theta]=meshgrid(phi,theta);
% plot(phi,theta,'b+',phip,thetap,'ko')
%
% plotcont([],[],3); hold on; lon=20; lat=30; 
% [X,Y,Z]=sph2cart(lon*pi/180,lat*pi/180,1);
% plot3(X,Y,Z,'o','MarkerF','b','MarkerE','b'); hold on;
% [thetap,phip,rotmats]=rottp(pi/2-lat*pi/180,lon*pi/180,...
%       10*pi/180,10*pi/180,0);
% [X,Y,Z]=sph2cart(phip,pi/2-thetap,1);
% plot3(X,Y,Z,'v','MarkerF','r','MarkerE','r'); hold off
% xlabel('x'); ylabel('y'); zlabel('z')
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/11/2010

% Supply defaults
defval('theta',linspace(0,pi,30))

if ~isstr(theta)
  defval('phi',linspace(0,2*pi,30))
  defval('alfa',rand(1)*2*pi)
  defval('bita',rand(1)*pi)
  defval('gama',rand(1)*2*pi)
  
  % Matrices unequal-size inputs 
  [ty,tx]=size(theta);
  [py,px]=size(phi);
  
  if ~all([tx ty]==[px py])
    [phi,theta]=meshgrid(phi(:),theta(:)');
    [ty,tx]=size(theta);
    [py,px]=size(phi);
  end

  % SPH2CART, see DT (A.106)
  x=sin(theta).*cos(phi);
  y=sin(theta).*sin(phi);
  z=cos(theta);

  % Successive PASSIVE rotations see DT (C.238)-(C.240) - vectorized
  rotmats=rotz(gama)*roty(bita)*rotz(alfa);

  % Note that this now is equivalent to this; keep checking in case I mess up
  difer(rotmats-rotz(gama+pi/2)*roty(pi/2)*rotz(0)...
	*rotz(bita)*roty(-pi/2)*rotz(alfa-pi/2),[],[],NaN);
  % This gives the coordinates of the original vector in the rotated system
  % so the intepretation is that it gives the coordinates of the
  % inversely-rotated vector in the original system... hence the
  % explanations, which are all backward from DT, though the operations
  % are all the same. (Note: flipping the signs is flipping the directions.)
  xyzp=rotmats*[x(:)'; y(:)' ; z(:)'];

  % CART2SPH, see DT (A.107)
  thetap=atan2(sqrt(xyzp(1,:).^2+xyzp(2,:).^2),xyzp(3,:));
  phip=atan2(xyzp(2,:),xyzp(1,:));

  % Return the same sizes as input
  thetap=reshape(thetap,[ty tx]);
  phip=reshape(phip,[py px]);
  
  % Return output if so desired
  varns={thetap,phip,rotmats};
  varargout=varns(1:nargout);
elseif strcmp(theta,'demo1')
  clf
  c11=[270 15];
  cmn=[330 -62];
  [a,b,lola]=plotcont(c11,cmn);
  % In degrees
  alp=10:5:60;
  bta=0;
  gam=0;
  for index=1:length(alp)
    [thetap,phip]=rottp(pi/2-lola(:,2)*pi/180,lola(:,1)*pi/180,...
			alp(index)*pi/180,bta*pi/180,gam*pi/180);
    hold on
    ch=plot([phip+(phip<0)*2*pi]*180/pi,90-thetap*180/pi,'b');
    title(sprintf('Points rotated over z by %i%s from y to x',...
		  alp(index),str2mat(176)))
    ylim([cmn(2) c11(2)])
    xlim([c11(1)-alp(index) cmn(1)])
    pause
    delete(ch)
  end
  hold off
elseif strcmp(theta,'demo2')
  clf
  c11=[270 15];
  cmn=[330 -62];
  [a,b,lola]=plotcont(c11,cmn);
  % In degrees
  alp=0;
  bta=10:5:60;
  gam=0;
  for index=1:length(bta)
    [thetap,phip]=rottp(pi/2-lola(:,2)*pi/180,lola(:,1)*pi/180,...
			alp*pi/180,bta(index)*pi/180,gam*pi/180);
    hold on
    ch=plot([phip+(phip<0)*2*pi]*180/pi,90-thetap*180/pi,'b');
    title(sprintf('Points rotated over y by %i%s from x to z',...
		  bta(index),str2mat(176)))
    ylim([cmn(2) c11(2)+bta(index)])
    xlim([c11(1) cmn(1)])
    pause
    delete(ch)
  end
  hold off
elseif strcmp(theta,'demo3')
  clf
  [a,b,lola]=plotcont;
  % In degrees
  alp=0;
  % We should start seeing Antarctica come in on the Greenwich i.e. the x
  % z plane and the North Pole come into center view
  bta=10:10:90;
  gam=0;
  for index=1:length(bta)
    [thetap,phip]=rottp(pi/2-lola(:,2)*pi/180,lola(:,1)*pi/180,...
			alp*pi/180,bta(index)*pi/180,gam*pi/180);
    hold on
    ch=plot([phip+(phip<0)*2*pi]*180/pi,90-thetap*180/pi,'b');
    title(sprintf('Points rotated over y by %i%s from x to z',...
		  bta(index),str2mat(176)))
%    ylim([cmn(2) c11(2)+bta(index)])
%    xlim([c11(1) cmn(1)])
    pause
    delete(ch)
  end
  hold off
end
