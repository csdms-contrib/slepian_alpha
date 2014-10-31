function [Xrot,Yrot]=rotdom(x,y,theta)
% [Xrot,Yrot]=ROTDOM(x,y,theta)
%
% Rotates the domain given by vectors 'x' and 'y' by an angle 'theta'
% in counterclockwise direction. 'x' and 'y' not out of meshgrid.
% 'theta' is in radians
%
% EXAMPLE:
%
% rotdom('demo');
%
% Written by FJS, August 6th 1998

if ~isstr(x)
  [X,Y]=meshgrid(x,y);
  [m,n]=size(X);
  rotmat=[cos(theta) -sin(theta) ; sin(theta) cos(theta)];
  
  XY=[X(:) Y(:)]*rotmat';
  
  Xrot=reshape(XY(:,1),m,n);
  Yrot=reshape(XY(:,2),m,n);
elseif strmatch(x,'demo')
  [X,Y]=meshgrid(-10:0.1:10); roto=rand*pi;
  [Xrot,Yrot]=rotdom(-10:0.1:10,-10:0.1:10,roto);
  sx=rand*10; sy=rand*10;
  subplot(211)
  pcolor(X,Y,gauss2(X,Y,sx,sy)) ; shading flat
  xlabel([ 'x ; \sigma_x= ',num2str(ceil(sx))]) 
  ylabel([ 'y ; \sigma_y= ',num2str(ceil(sy))]) 
  title([ 'Original function']) ; axis ij ; axis image
  subplot(212)
  pcolor(X,Y,gauss2(Xrot,Yrot,sx,sy)) ; shading flat 
  xlabel('x') ; ylabel('y') ; axis ij ; axis image
  title([ 'Rotated by ',num2str(ceil(roto*180/pi)),'\circ ccw'])
end


