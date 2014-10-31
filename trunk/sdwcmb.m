function sdwcmb
% SDWCMB
%
% Simons, Dahlen and Wieczorek (2005)
% Polar and antipolar caps 
%
% Last modified by fjsimons-at-alum.mit.edu, August 19th, 2004

TH=40;
m=0;
nth=32;
L=18;
np=4;

% T1 is the top cap
[E,V,th,C,T]=grunbaum(TH,L,m,nth);
[E2,V2,jk1,jk2,jk3,jk4,jk5,jk6,jk7,jk8,K]=sdwcap(TH,L,m,nth,-1);

% Grunbaum's nominal eigenvalues
VG=diag(C'*K*C);

[ah,ha]=krijetem(subnum(3,4));

% Now try to add a 3-D collar around this thing
[x1,y1,z1]=latitude(90-TH,1);
[x2,y2,z2]=latitude(90-TH,1);

% This doesn't work, really
%for index=1:size(E,2)
%  F=scale(E(:,index),[0 1]);
%  E(abs(F)<1/100,index)=NaN;
%end

for index=1:4
  axes(ah(index))
  plotonsphere(E(:,index)*ones(1,nth),0.15)  
  view(0,30)
  tits=sprintf('%s_%i = %8.3e','\lambda',index,VG(index));
  t(index)=title([sprintf('%s%s^{%s}',pref(tits,'e'),' \times 10',suf(tits,'e'))]);
  hold on
%  plot3(x1,y1,z1,'w')
  EM(index)=max(abs(E(:,index)));
end

for index=5:8
  axes(ah(index))
  plotonsphere(-E(:,L-7+index)*ones(1,nth),0.15)  
  view(0,-30)
  tits=sprintf('%s_{%i} = %8.3e','\lambda',L-7+index,VG(L-7+index));
  t(index)=title([sprintf('%s%s^{%s}',pref(tits,'e'),' \times 10',suf(tits,'e'))]);
  hold on
%  plot3(x2,y2,z2,'w')
  EM(index)=max(abs(E(:,index)));
end

for index=9:12
  axes(ah(index))
  plotonsphere(-E2(:,L-11+index)*ones(1,nth),0.15)  
  view(0,-30)
  tits=sprintf('%s_{%i} = %8.3e','\lambda',L-11+index,V2(L-11+index));
  t(index)=title([sprintf('%s%s^{%s}',pref(tits,'e'),' \times 10',suf(tits,'e'))]);
  hold on
%  plot3(x2,y2,z2,'w')
  EM(index)=max(abs(E(:,index)));
end

kelicol
for index=1:12
  axes(ah(index))
  shading faceted  
  caxis([-EM(index) +EM(index)])
end

movev(t(1:4),-0.5)
movev(t(5:12),0.5)

delete(ah(9:12))

fig2print(gcf,'portrait')
figdisp([],[],'-painters')

