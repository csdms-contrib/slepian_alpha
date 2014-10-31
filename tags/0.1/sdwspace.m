function sdwspace
% SDWSPACE
%
% Simons, Dahlen and Wieczorek (2006) - Figure 5.1
% Plots spatial functions for spherical polar cap
% for various angular orders with eigenvalue grey shading. 
%
% Last modified by fjsimons-at-alum.mit.edu, August 19th, 2004

TH=40; TH=20
nth=128;
L=18; SN=4;
M=4;

clf
[ah,ha]=krijetem(subnum(M+1,SN));

for m=0:M
  [E1,V1{m+1},N1,th1,C1,ngl1,ngl2]=sdwcap(TH,L,m,nth);
  % Note that we are requesting but not using C2, never mind the error 
  [E2,V2{m+1},N2,th2,C2]=sdwcapt(TH,L,m,nth);
  E1m(m+1)=max(abs(E1(:)));
  E2m(m+1)=max(abs(E2(:)));
  for ondex=1:SN
    axes(ah(m*SN+ondex))
    % pm(m+1,ondex)=plot(th1,E1(:,ondex),'k-');
    % hold on
    % pg(m+1,ondex)=plot(th2,E2(:,ondex),'-','Color',grey(6));
    pm(m+1,ondex)=plot(th1,E1(:,ondex),'-','Color',grey(6));
    hold on
    pg(m+1,ondex)=plot(th2,E1(1:size(E2,1),ondex),'-','Color','k');
    drawnow
    % set(ah(m*SN+ondex),'Color',grey(V1(ondex)*10))
    % if V1(ondex)<0.4
    %  set(pm(m+1,ondex),'Color','w')
    % end
  end
end

set(ah,'xlim',[0 180],'xtick',[0 TH 180],'xgrid','on','ygrid','on',...
       'ylim',[-1.1 1.1]*max([E1m E2m]))
nolabels(ah(1:end-SN),1)
nolabels(ha(M+2:end),2)
%nolabels(ha(1:M+1),2)
longticks(ah,1/2)

% Eigenvalue labels
nf=9;
for m=0:M
  for ondex=1:SN
    axes(ah(m*SN+ondex))
        t{m+1,ondex}=sprintf('%s = %9.6f','\lambda',V1{m+1}(ondex));
	if V1{m+1}(ondex)<0.975
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
  deggies(ah(index),1)
  xl1(index)=xlabel(sprintf('colatitude %s','\theta'));
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

set([ah],'Fonts',nf)

fig2print(gcf,'portrait')

% set(gcf,'Color','w','Inv','off')
shrink(ah,1,1/1.1)

figdisp 
