function sdwfried
% SDWFRIED
%
% Makes fried-egg plots for +- angular orders
% Figure 5.4 from  Simons et al., SIAM Review (2006)
%
% Run twice. Second time round the color axis is OK
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2009

TH=40;
L=18;

[lrnk,mrnk,lval]=sdwelm(TH,L);

% Nearly double the amount of requested tapers
lrnk=gamini(lrnk,(mrnk~=0)+1);
lval=gamini(lval,(mrnk~=0)+1);
mrnk=gamini(mrnk,(mrnk~=0)+1);
% This trick only works for the single cap
mrnk((diff(mrnk)==0))=-mrnk((diff(mrnk)==0));

[ah,ha]=krijetem(subnum(8,4));
for index=1:length(ah) 
  fnpl=sprintf('%s/SDWFRIED-%i-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'LOCALIZE'),TH,L,index);
  if exist(fnpl,'file')==2
    load(fnpl)
    disp(sprintf('%s loaded by SDWFRIED',fnpl))
    axes(ah(index))
    cbb=gray(10);
    cbb=kelicol;
    ispl=data/max(abs(data(:)));
    ispl(abs(ispl)<1/100)=NaN;
    imagefnan([-1 1],[1 -1],ispl,cbb,[-1 1])
    colormap(cbb); hold on
    ax2(1)=circ(1); ax2(2)=circ(sin(TH*pi/180));
    set(ax2(1:2),'LineW',1)
    set(ax2(2),'LineS','--')
    axis([-1.0100    1.0100   -1.0100    1.0100])
    axis off
    if index==0
      [axlim,handl,XYZ]=plotcont([0 90],[360 -90],4);
      delete(handl)
      hold on
      p1=plot(XYZ(:,1),XYZ(:,2),'k-','LineWidth',1);
    end
    tl(index)=title(sprintf('%s_%i = %6.3f ; m = %i',...
			    '\lambda',...
			    lrnk(index),lval(index),mrnk(index)));
    drawnow
  else
    [E,V,N,th,C]=sdwcap(TH,L,mrnk(index),128,[],2);
    E=E{lrnk(index)};
    Emax=max(E(:)); thresh=100;
    Ed=E;
    Ed(abs(Ed)<Emax/thresh)=0;
    axes(ah(index))
    [data,ch,ph]=plotplm(Ed,[],[],5,[],TH);
    tl(index)=...
	title(sprintf('%s = %6.3f ; m = %i','\lambda',...
		      lval(index),mrnk(index)));
    delete(ch)
    drawnow
    save(fnpl,'data')
  end
end
seemax(ah,3)
fig2print(gcf,'tall')
set(ah,'CameraV',6.25)
% Last-minute cosmetics
movev(tl,-0.25)

figdisp


