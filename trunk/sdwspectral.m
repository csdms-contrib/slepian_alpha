function sdwspectral 
% SDWSPECTRAL
%
% Simons, Dahlen and Wieczorek (2006) - Figure 5.2
% Plots spectral functions for spherical polar cap
% for various angular orders.
%
% Last modified by fjsimons-at-alum.mit.edu, 04/06/2009

% Check absolute power; then don't do deciBel since the maxima may not be
% reached  
TH=40;
nth=128;
L=18; SN=4;
M=4;
Lnyq=nth-1;

clf
[ah,ha]=krijetem(subnum(M+1,SN));
yls=[-140 10];

for m=0:M
  [E1,V1{m+1},N1,th1,C1,ngl1,ngl2]=sdwcap(TH,L,m,nth);
  [E2,V2{m+1},N2,th2,C2]=sdwcapt(TH,L,m,nth);  
  
  E1m(m+1)=max(abs(E1(:)));
  E2m(m+1)=max(abs(E2(:)));
  narm2=1./(2*[m:Lnyq]'+1);
  narm2=1;
  narm1=1./(2*[m:L]'+1);
  narm1=1;

  for ondex=1:SN
    axes(ah(m*SN+ondex))
    spec=decibel(narm2.*C2(:,ondex).^2);
    pg(m+1,ondex)=plot([m:Lnyq],spec,'-','Color',grey(6));
    hold on
    pm(m+1,ondex)=plot([m:L],spec(1:L-m+1),'-','Color','k');
    drawnow
    plot([m m],yls,'k:')
    hold off
  end
end

set(ah,'xlim',[0 127],'xtick',[0 L Lnyq],...
       'xgrid','on','ygrid','on','ylim',yls)
nolabels(ah(1:end-SN),1)
nolabels(ha(M+2:end),2)
longticks(ah,1/2)

% Eigenvalue labels
nf=9;
for m=0:M
  for ondex=1:SN
    axes(ah(m*SN+ondex))
        t{m+1,ondex}=sprintf('%s = %9.6f','\lambda',V2{m+1}(ondex));
	if V2{m+1}(ondex)<0.975
	  lox='lr';
	else
	  lox='ur';
	end
	[bh(m+1,ondex),th(m+1,ondex)]=...
	    boxtex(lox,ah(m*SN+ondex),t{m+1,ondex},nf);
  end  
end
set(th,'FontS',nf-1)

for index=SN*M+1:SN*(M+1)
  axes(ah(index))
  xl1(index)=xlabel(sprintf('degree l'));
end

for index=1:M+1
  serre(ah([1:SN]+(index-1)*SN),1/3)
  axes(ha(index))
  ylb(index)=ylabel(sprintf('m = %i',index-1));
end

for index=1:SN
  axes(ah(index))
  tlb(index)=title(sprintf('%s = %i','\alpha',index));
end

set([pm(:) ; pg(:)],'LineW',1)

set(ah,'Fonts',nf)

axes(ha(M+1))
%tt=text(-60,-195,'dB','FontS',nf);
tt=text(-60,-170,'dB','FontS',nf);

fig2print(gcf,'portrait')

% set(gcf,'Color','w','Inv','off')
shrink(ah,1,1/1.1)

figdisp 
