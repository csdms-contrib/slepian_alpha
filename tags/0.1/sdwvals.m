function sdwvals
% SDWVALS
%
% Eigenvalue structure for various concentration regions.
% Mixed-order ranking.
%
% Last modified by fjsimons-at-alum.mit.edu, 04/24/2009

TH=[10 20 30 40];
L=18;

clf
[ah,ha,H]=krijetem(subnum(2,2));

yls=[-0.1 1.1];

symbs={'o','x','s','+','v','*','^','d','<','>','p','h',...
       'o','x','s','+','v','*','^','d'};
xmax=60;

more off
for index=1:length(TH)
  legsi{index}=sprintf('%s = %i%s ','\Theta',TH(index),str2mat(176));
  SN(index)=floor(L/180*TH(index));
  Nall(index)=(L+1)^2*(1-cos(TH(index)*pi/180))/2;

  [lrnk,mrnk,lval,VV,Vsum]=sdwelm(TH(index),L);
  
  ldubs=gamini(lval,~~mrnk+1);
  mrnk=gamini(mrnk,~~mrnk+1);
  axes(ah(index))
  for ondi=1:length(ldubs)
    p(ondi,index)=plot(ondi,ldubs(ondi),symbs{mrnk(ondi)+1});
    hold on
  end
  hold on
  plot(round([Nall(index) Nall(index)]),yls,'k:')
  plot([0 xmax],[0.5 0.5],'k:')
  plot([0 xmax],[0 0],'k:')
  plot([0 xmax],[1 1],'k:')
  set(ah(index),'xlim',[0 xmax],'ylim',yls,'xgrid','off','ygrid','off',...
		'xtickl',[1 10:10:xmax],'xtick',[1 10:10:xmax],...
		'ytick',[0:0.25:1])
  [bh(index),th(index)]=boxtex('ur',ah(index),legsi{index},12);
  drawnow
end

axes(ah(4))
fb=fillbox([2 18 0.88 -0.05],'w');
for ondi=1:12
  ypo=0+0.075*(ondi-1);
  pl(ondi,1)=plot(4,ypo,symbs{ondi});
  hold on
  tl(ondi,1)=text(7,ypo,sprintf('m = %s %i','\pm',ondi-1),'FontS',8);
end

% Now make the plot beautiful
movev(th,-.00)
longticks(ah)
set([p(~~p(:)) ; pl(~~pl(:))],'MarkerS',4,'MarkerF',grey,'MarkerE','k')
axes(ha(1))
al(1)=ylabel('eigenvalue \lambda');
axes(ha(2))
al(2)=ylabel('eigenvalue \lambda');
axes(ha(2))
xl(1)=xlabel('rank');
axes(ha(4))
xl(2)=xlabel('rank');

nolabels(ha(3:4),2)
nolabels(ah(1:2),1)

serre(H',1/2,'down')
serre(H,1/2,'across')

for ind=1:4
  xx(ind)=xtraxis(ah(ind),round(Nall(ind)),...
		  {sprintf('N = %i',round(Nall(ind)))});
end
longticks(xx)

set([xl al],'FontS',13)
set([ ah],'FontS',12)

figdisp


