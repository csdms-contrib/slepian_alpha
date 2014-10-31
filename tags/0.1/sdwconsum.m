function sdwconsum(L)
% SDWCONSUM(L)
%
% Plots the eigenvalue-weighted sum of all Slepian functions localized to
% all of the continents except Antarctica, as in: 
% Simons, Dahlen and Wieczorek (SIAM, 2006), Figure 6.4
%
% INPUT:
%
% L       Spherical harmonic degree of the bandwidth [default: 18]
%
% Last modified by fjsimons-at-alum.mit.edu, 03/25/2009

doms={'africa', 'eurasia', 'namerica', 'australia', 'greenland', ...
      'samerica'};
sa=0;
for index=1:length(doms)
  sa=sa+spharea(doms{index});
end
defval('L',18);

% Calculate the Shannon number
N=(L+1)^2*sa;

% Specify the truncation levels in the sum
ens=ceil([N/4 N/2 N (L+1)^2]);
legsi{1}=sprintf('1%s N/4','\rightarrow');
legsi{2}=sprintf('1%s N/2','\rightarrow');
legsi{3}=sprintf('1%s N','\rightarrow');
legsi{4}=sprintf('1%s (L+1)^2','\rightarrow');

clf
[ah,ha]=krijetem(subnum(2,2));

for index=1:length(ens)
  axes(ah(index))
  fnpl=sprintf('%s/SDWALLCONS-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'KERNELC'),L,ens(index));
  if exist(fnpl,'file')==2
    load(fnpl)
  else
    F=sdwallcons(doms,L,0);
  end
  F=scale(F,[0 1]);
  F(F<max(F(:))/100)=NaN;
  % imagefnan([0 90],[360 -90],F,'default',minmax(F))
  % imagefnan([0 90],[360 -90],F,'kelicol',minmax(F))
  imagefnan([0 90],[360 -90],F,'kelicol',[-1 1])
  set(gca,'ytick',[-90:45:90])
  set(gca,'xtick',[0:90:360])
  deggies(gca)
  [jk,a]=plotcont; set(a,'linew',1)
  axis image
  [bh(index),th(index)]=boxtex('ll',ah(index),legsi{index},14);
end

% Cosmetics
longticks(ah,3/2)
nolabels(ha(3:4),2)
nolabels(ah(1:2),1)
set(ah,'Camerav',6.5)
serre(ha(1:2),1.25,'down')
serre(ha(3:4),1.25,'down')

% TRICK TO GET A COLOR BAR
%a=colorbarf('hor',10,'Helv',[0.25 0.05 0.5 0.03]);
caxcon=[0 1];
caxoc=[-1 0];
h=axes;
[cb,xcb]=addcb('hor',caxcon,caxoc,'kelicol');
delete(h)
shrink(cb,2,2)
movev(cb,-0.075)
set(xcb,'String','cumulative energy')
axes(cb)
xlim([0 1])

fs=15;
set(cb,'xtick',[0 1],'xtickl',{'0' 'N/A'},'FontS',fs)
set(xcb,'FontS',fs)
set(ah,'FontS',fs-1)

movev([ah cb],.1)

fig2print(gcf,'landscape')
figdisp
