function sdwsumall
% SDWSUMALL
%
% Simons, Dahlen and Wieczorek (2006)
% Sums all eigenfunctions.
%
% First calculate a bunch of cumulative sums. See SDSUMALL and SDSUMALL2
%
% Should still modify; we were originally plotting more than one
% longitude. See SDSUMALL where it is verified.
%
% Last modified by fjsimons-at-alum.mit.edu, 09.08.2005
% Stupid program, does double work. Later.

TH=[10 20 30 40];
L=18;
nth=128;

for ondex=1:length(TH)
  % Figure out order to sum them in 
  [lrnk,mrnk,lval]=sdwelm(TH(ondex),L);

  % Nearly double the amount of requested tapers for repeated abs(m)
  lval=gamini(lval,(mrnk~=0)+1);
  lrnk=gamini(lrnk,(mrnk~=0)+1);
  mrnk=gamini(mrnk,(mrnk~=0)+1);
  mrnk((diff(mrnk)==0))=-mrnk((diff(mrnk)==0));

  % This condition no good here, since SDWCAP has VCUT on
  %  if length(lrnk)~=(L+1)^2
  %    error('Something wrong')
  %  end
  
  % Calculate the Shannon number
  N(ondex)=round((L+1)^2*(1-cos(TH(ondex)*pi/180))/2);

  disp(sprintf('TH = %i ; Shannon number is %i',TH(ondex),N(ondex)))
  
  % Must have all the pairs of orders to be longitudinally independent
  if abs(mrnk(N(ondex)))==abs(mrnk(N(ondex)+1))
    N(ondex)=N(ondex)+1;
  end
  
  % Calculate cumulative sums weighted by eigenvalue!
  F=0; 
  for index=1:length(lrnk)
    % Note this is not very efficient since all the nonzero orders are
    % simpy done twice!
    fnpl=sprintf('%s/SDWSUMCUM-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'LOCALIZE'),...
		 TH(ondex),L,index);
    if exist(fnpl,'file')~=2
      [E,V,Ns,th,C]=sdwcap(TH(ondex),L,mrnk(index),nth,[],2);
      % We are summing the FULL longitudinal matrix, i.e. with the
      % sqrt(2)cos(mphi) and sqrt(2)sin(mphi) already in there!
      F=F+(E{lrnk(index)}).^2.*V(lrnk(index));
      eval(sprintf('save %s F',fnpl))
    end
  end
end

% Now plot them
clf
[ah,ha]=krijetem(subnum(2,2));
xls=[0 180];
yls=[-1 35];
cols=[grey(3) ; grey(6) ; 0 0 0 ];
NA=(L+1)^2/4/pi;
for ondex=1:length(TH)
  axes(ah(ondex))
  legsi{ondex}=sprintf(' %s = %i%s  ','\Theta',TH(ondex),str2mat(176)); 
  [lrnk,mrnk,lval]=sdwelm(TH(ondex),L);
  lrnk=gamini(lrnk,(mrnk~=0)+1);
  pott=[length(lrnk) N(ondex)];
  post{ondex}=...
     {'1\rightarrow N' '1\rightarrow (L+1)^2' 'unweightedi'};
  for index=1:length(pott)
    fnpl=sprintf('%s/SDWSUMCUM-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'LOCALIZE'),...
		 TH(ondex),L,pott(index));
    eval(sprintf('load %s',fnpl))
    noteight=8;
    p{index}(:,ondex)=plot(linspace(0,180,size(F,1)),...
	 F(:,round(linspace(1,size(F,2),noteight))),...
	 'Color',cols(index,:));
    hold on
  end  
  set(ah(ondex),'xtick',[0 TH(ondex) 180],'xgrid','on',...
		'ytick',[0 NA],'ygrid','off')
  set(ah,'xlim',xls,'ylim',yls)
  % ADD THE UNWEIGHTED STUFF
  G=0;
  % Sum over all m's WITHOUT THE VCUT! But it still reads the saved ones
  % first, so this is also no good, really
  for jndex=0:L
    [E2,V2,N2,th]=sdwcap(TH(ondex),L,jndex,nth,-1);
    % Remember the factor of sqrt(2-dom) was already in there!
    G=G+sum(E2.^2,2);
  end
  % Check this outcome
  if any(abs(G-NA)>1e-10)
    warning(sprintf(...
	'Unweighted sum exceeds N/A for L= %i and TH= %i',...
	L,TH(ondex)))
  end

  punw(ondex)=plot(th,G,'LineW',1,'Color',cols(3,:),'LineS','--');
  % END ADD UNWEIGHTED STUFF
  [bhh(ondex),thh(ondex)]=boxtex('ur',ah(ondex),legsi{ondex},12);
  drawnow
end
longticks(ah)

axes(ha(1))
yl(1)=ylabel('cumulative energy');
axes(ha(2))
yl(2)=ylabel('cumulative energy');
axes(ah(3))
xl(1)=xlabel('colatitude \theta');
axes(ah(4))
xl(2)=xlabel('colatitude \theta');
deggies(ah(2:4),1)
set(ah(1),'xtick',[TH(1) 180],'xtickl',[TH(1) 180])
deggies(ah(1),1)
% This needs to be after deggies
set(ha(1:2),'ytickl',{'' ''})
set(ha(3:4),'ytickl',{'0' 'N/A'},'yaxisl','right')
set([yl xl],'FontS',13)
set(ah,'FontS',12)

set([p{:}],'LineW',1)
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')

for index=1
  axes(ah(index))
  loh(index)=legend(post{index});  
  set(getkids(loh(index),2),'Color',cols(3,:),'LineS','--')
  set(getkids(loh(index),5),'Color',cols(1,:))
  set(getkids(loh(index),8),'Color',cols(2,:))
end
movev(loh,-.070)
moveh(thh(:),5)

set(findobj('string','unweightedi'),'string','unweighted')

fig2print(gcf,'portrait')
figdisp

