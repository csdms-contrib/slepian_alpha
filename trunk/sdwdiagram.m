function sdwdiagram
% SDWDIAGRAM
%
% Makes a diagram of the spherical set-up of the program
% Simons, Dahlen and Wieczorek, Figure 1.
%
% Last modified by fjsimons-at-alum.mit.edu, August 19th, 2004

% Which vector to plot
ang=40;
% Down to this z level for the projection
lz=-0.2;
% Rotation of the geodesic
rotg=-30;
% Rotation of the X-axis
rots=10;

clf
[ah,ha]=krijetem(subnum(2,2));

% FIRST FIGURE PANEL ------------------------------
axes(ah(1))
[h,cord]=circ(1,[-pi/2 pi/2]); delete(h)
eq(1)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
hold on
[h,cord]=circ(1,[0 3*pi/2]); delete(h)
eqd(2)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),':');
[h,cord]=circ(1); delete(h)
ol(1)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));

% Plot the X-axis
[axt(1),axh(1)]=plotx;

% Plot the Y-axis
[axt(2),axh(2)]=ploty;

% Plot the Z-axis
[axt(3),axh(3)]=plotz;

% Plot the random vector
ft=[1 1]; 
vax=arrow(0,0,0,0.85,'v',2,ang,ft);
xdv=get(vax,'Xdata'); ydv=get(vax,'Ydata'); delete(vax)
axt(4)=plot3(zeros(size(ydv{1})),xdv{1}-ft(1),ydv{1}-ft(2));
%axh(4)=plot3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2));
ydv{2}(4)=ydv{2}(1);
xdv{2}(4)=xdv{2}(1);
axh(4)=fill3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2),'k');
% Plot first arclength
[h,cord]=circ(0.5,[pi/2-ang/180*pi pi/2]); delete(h)
arcl(1)=plot3(zeros(size(cord(:,1))),cord(:,1),cord(:,2));

% Plot projections
proj(1)=plot3([0 0],repmat(max(xdv{1})-ft(1),1,2),...
	      [max(ydv{1})-ft(1) lz]);
proj(2)=plot3([0 0],[0 max(xdv{1})-ft(1)],[0 lz]);

% Plot second arclength
[h,cord]=circ(0.4,[-rots*pi/180 ang/180*pi]); delete(h)
arcl(2)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
onax=axis; onax=onax.*1.1;
viewpars(onax)

% SECOND FIGURE PANEL ------------------------------
axes(ah(2))
% Plot Earth outline
[h,cord]=circ(1); delete(h)
ol(2)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
hold on
% Plot tilted equatorial plane
[h,cord]=circ(0.985,[-pi/2 pi/2]); delete(h)
rot2=-20;
cro=rotx(rot2*pi/180);
cord=[cro*[cord zeros(size(cord(:,1)))]']';
eq(2)=plot3(cord(:,1),cord(:,2),cord(:,3));
[h,cord]=circ(0.975,[0 3*pi/2]); delete(h)
cro=rotx(rot2*pi/180);
cord=[cro*[cord zeros(size(cord(:,1)))]']';
eqd(3)=plot3(cord(:,1),cord(:,2),cord(:,3),':');
% Plot two vectors along this great circle
ang=70;
%vax=arrow(0,0,0,0.85,'v',2,ang,ft);
vax=arrow(0,0,0,1.1,'v',2,ang,ft);
xdv=get(vax,'Xdata'); ydv=get(vax,'Ydata'); delete(vax)
axt(8)=plot3(zeros(size(ydv{1})),xdv{1}-ft(1),ydv{1}-ft(2));
%axh(8)=plot3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2));
ydv{2}(4)=ydv{2}(1);
xdv{2}(4)=xdv{2}(1);
axh(8)=fill3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2),'k');
ang=200;
vax=arrow(0,0,0,0.44,'v',1,ang,ft);
xdv=get(vax,'Xdata'); ydv=get(vax,'Ydata'); delete(vax)
axt(9)=plot3(zeros(size(ydv{1})),xdv{1}-ft(1),ydv{1}-ft(2));
%axh(9)=plot3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2));
ydv{2}(4)=ydv{2}(1);
xdv{2}(4)=xdv{2}(1);
axh(9)=fill3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2),'k');
% Plot arclength connecting two vectors
[h,cord]=circ(0.4,[-7.5 87.5]*pi/180); delete(h)
cro=rotx(rot2*pi/180);
cord=[cro*[cord zeros(size(cord(:,1)))]']';
arcl(4)=plot3(cord(:,1),cord(:,2),cord(:,3));
viewpars(onax)

% THIRD FIGURE PANEL
axes(ah(3))
% Plot Earth outline
[h,cord]=circ(1); delete(h)
ol(3)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
hold on
[h,cord]=circ(1,[-pi/2 pi/2]); delete(h)
eq(5)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
[h,cord]=circ(1,[0 3*pi/2]); delete(h)
eqd(4)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),':');
% First blob
[lo,la]=blob; 
lo=scale(lo,[0 60]*pi/180);
la=scale(la,[20 50]*pi/180);
[X,Y,Z]=sph2cart(lo,la,1);
arcl(3)=fill3(X,Y,Z,grey);
% Second blob
[lo,la]=blob; 
lo=scale(lo,[-30 0]*pi/180);
la=scale(la,[10 -30]*pi/180);
[X,Y,Z]=sph2cart(lo,la,1);
arcl(6)=fill3(X,Y,Z,grey);
curs=0;
if curs==1
  [phint,thp,php]=phicurve([pi/2-la(:) lo(:)],...
			   linspace(pi/2-max(la),pi/2-min(la),10));
  % Now get the great circle coordinates
  thp=thp'; php=php';
  for index=1:size(php,1)
    lola=grcircle(...
	[php(index,1) pi/2-thp(index,1)],...
	[php(index,2) pi/2-thp(index,2)],100);
    [X,Y,Z]=sph2cart(lola(:,1),lola(:,2),1);
    hasj(index)=plot3(X,Y,Z,'k');
  end
end
viewpars(onax)

% FOURTH FIGURE PANEL
axes(ah(4))
% Plot Earth outline
[h,cord]=circ(1); delete(h)
ol(4)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
hold on
viewpars(onax)
[axt(5),axh(5)]=plotx;
[axt(6),axh(6)]=ploty;
[axt(7),axh(7)]=plotz;
[h,cord]=circ(1,[-pi/2 pi/2]); delete(h)
eq(3)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
[h,cord]=circ(1,[0 3*pi/2]); delete(h)
eqd(1)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),':');
% Plot polar cap
ang=60;
[h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
eq(4)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
[h,cord2]=circ(1,[ang 180-ang]*pi/180); delete(h)
eqx=fill3([cord(:,1) ;  ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ;  ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang*pi/180) ; ... 
	     ; cord2(:,2)]',...
	    grey);

% Indicate the angle
proj(3)=plot3([0 0],[0 cos(ang*pi/180)],[0 sin(ang*pi/180)]);
proj(4)=plot3([0 0],[0 -cos(ang*pi/180)],[0 sin(ang*pi/180)]);
[h,cord]=circ(0.5,[pi/2-ang/180*pi/2 pi/2]); delete(h)
arcl(5)=plot3(zeros(size(cord(:,1))),cord(:,1),cord(:,2));
if curs==1
  % Plot hashings
  ang=65;
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  eqh(1)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  ang=70;
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  eqh(2)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  ang=75;
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  eqh(3)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  ang=80;
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  eqh(4)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  ang=85;
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  eqh(5)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
else
  eqh=[];
end
viewpars(onax)

% Cosmetics
serre(ah(1:2),1,'across')
serre(ah(3:4),1,'across')
serre(ha(1:2),0.6,'down')
serre(ha(3:4),0.6,'down')
set([eq eqh eqd ol axt proj],'LineW',1,'Color','k')
set([arcl([1 2 4 5])],'LineW',1,'Color','k')
set(eqh,'LineW',0.5)
set([eq proj],'LineS',':')
set(ah,'camerav',5.5) 

% Plot labels
axes(ah(1))
lpost=[1.9  -0.4    0  ;
      0.5    1.2    0  ;
      0      0.175  1.2;
      -0.2   1.05   0.3;
      0.4   -0.2   -0.4;
      0      0.5    0.7];
axpo=[1 1 1 2 2 1];
ltxt={'\bfx','\bfy','\bfz','\bfr ''','\bfr','\bfr'};
for i=1:6
  axes(ah(axpo(i)))
  l(i)=text(lpost(i,1),lpost(i,2),lpost(i,3),ltxt{i});
  % Extent only works in 2D
  % a=text(0,0,'x','fonts',get(l(i),'fontsize')); 
  % xte=get(a,'Extent');
  % delete(a)
  % Put the hat on the xi

  th(i)=text(lpost(i,1),lpost(i,2),lpost(i,3)+0.05,'\bf\^');
end
%set(th,'Horizontala','center')
axes(ah(1))
l(8)=text(0,0.2,0.6,'\theta');
l(9)=text(0.65,0.2,0.02,'\phi');
l(7)=text(0,0.4,-0.7,'\Omega');
axes(ah(2))
l(10)=text(0.3,0.4,0.125,'\Delta');
axes(ah(3))
l(11)=text(0,0.7,-0.05,'R_1');
l(13)=text(0,0.1,-0.6,'R_2');
axes(ah(4))
l(12)=text(0,0.125,0.575,'\Theta');

fig2print(gcf,'portrait')
figdisp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewpars(onix)
defval('onix',[-2 2 -2 2])
% Set viewing parameters
% Really would need another rotation around x
view(90,17.5); axis equal
axis(onix); axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h]=plotx
% Plots x-axis for the unit sphere
ft=[1 1]; % To get it on the right plane
% Rotation of the X-axis
rots=10;
xax=arrow(0,0,1.75,0,'h',2,rots,ft);
xdx=get(xax,'Xdata'); ydx=get(xax,'Ydata'); delete(xax)
t(1)=plot3(xdx{1}-ft(1),ydx{1}-ft(2),zeros(size(ydx{1})));
% Regular arrow head
% h(1)=plot3(xdx{2}-ft(1),ydx{2}-ft(2),zeros(size(ydx{2})));
% Filled arrow head, but unrotated
ydx{2}(4)=ydx{2}(1);
xdx{2}(4)=xdx{2}(1);
h(1)=fill3(xdx{2}-ft(1),ydx{2}-ft(2),zeros(size(ydx{2})),'k');
% Filled and properly rotated arrow head
% rots1=rotz(-rots*pi/180);
% rots2=rotx(pi/2);
% rots3=rotz(rots*pi/180);
% done=(rots3*rots2*rots1*[xdx{2}-ft(1),ydx{2}-ft(2),zeros(size(ydx{2}))]')';
% h(1)=fill3(done(:,1),done(:,2),done(:,3),'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h]=ploty
% Plots y-axis for the unit sphere
yax=arrow(0,0,0,1.4,'v',3);
xdy=get(yax,'Xdata'); ydy=get(yax,'Ydata'); delete(yax)
t(1)=plot3(xdy{1},ydy{1},zeros(size(ydy{1})));
% h(1)=plot3(zeros(size(ydy{2})),ydy{2},xdy{2});
ydy{2}(4)=ydy{2}(1);
xdy{2}(4)=xdy{2}(1);
h(1)=fill3(zeros(size(ydy{2})),ydy{2},xdy{2},'k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h]=plotz
% Plots z-axis for the unit sphere
zax=arrow(0,0,0,1.4,'v',3);
xdz=get(zax,'Xdata'); ydz=get(zax,'Ydata'); delete(zax)
t(1)=plot3(zeros(size(ydz{1})),xdz{1},ydz{1});
% h(1)=plot3(zeros(size(ydz{2})),xdz{2},ydz{2});
ydz{2}(4)=ydz{2}(1);
xdz{2}(4)=xdz{2}(1);
h(1)=fill3(zeros(size(ydz{2})),xdz{2},ydz{2},'k');

