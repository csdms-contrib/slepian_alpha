function varargout=sdwregions(par,fi,L)
% [ah,ha,bh,th]=SDWREGIONS(par,fi,L)
%
% INPUT:
%
% par     1  Eigenvalue plot [Simons et al. SIAM Review (2006) Fig 6.1]
%         3  Australia [Simons et al. SIAM Review (2006) Fig 6.2]
%         4  North America [Simons & Dahlen SPIE (2007) Fig 2]
%         5  Africa [Simons et al. SIAM Review (2006) Fig 6.3]
%         6  Eurasia
%         7  South America [Simons & Dahlen SPIE (2007) Fig 3]
%         8  Amazon basin
%         9  Orinoco basin
%        10  Greenland
%        11  "GPS North America", some arbitrary region
%        12  Antarctica [Simons Handbook of Geomathematics (2010) Fig 3]
%        13  The world's oceans [Slobbe et al. J Geodesy (2011) Fig 2]
% fi      0  On 4x3
%         1  On 3x4 
% L      Bandwidth (maximal spherical harmonic degree)
%
% See also KERNELC, LOCALIZATION, SDWALLCONS
%
% Last modified by fjsimons-at-alum.mit.edu, 09/25/2011

defval('par',1)
defval('fi',0)
defval('L',36)

clf
switch par
 case 1
  % These only for the eigenvalue plots; others have to add explicitly
  regs={'greenland','australia','namerica','africa','eurasia'};
  legs={'Greenland','Australia','N America','Africa','Eurasia'};

  L=[6 12 18 24]; % For the SIAM Review figure
  colt=linspace(0,10,length(regs));
  cols=repmat(0,1,length(regs));

  [ah,ha]=krijetem(subnum(2,2));
  
  for ondex=1:length(L)
    clear V C N
    for index=1:length(regs)
      N(index)=(L(ondex)+1)^2*spharea(regs{index});
      [V(:,index),C]=localization(L(ondex),regs{index});
    end
    axes(ah(ondex))
    for index=1:length(regs)
      pv(ondex,index)=plot(V(:,index),'-','Marker',symbol(index,1));
      hold on
      colo=grey(cols(index));
      colb=grey(colt(index));
      % set(pv(ondex,index),'Color',colo,'MarkerF',colo,'MarkerE',colo);
      set(pv(ondex,index),'Color',colo,'MarkerF',colb,'MarkerE','k');
    end
    xhi=1.5*max(N);
    plot([0 xhi],[0.5 0.5],'k:')
    plot([0 xhi],[0 0],'k:')
    plot([0 xhi],[1 1],'k:')

    if ondex==1
      xhi=13;
    end
    set(pv,'LineW',1)
    set(ah(ondex),'xtick',round(sort(N)),'xtickl',round(sort(N)),...
		  'xgrid','on','xlim',[0 xhi],'ytick',[0:0.25:1])
    tt(ondex)=text(xhi,-0.125,sprintf('%s %i','\rightarrow',...
		   (L(ondex)+1)^2),'Horiz','right');
    drawnow
  end
  
  axes(ah(1))
  yl(1)=ylabel('eigenvalue \lambda');
  axes(ah(3))
  xl(1)=xlabel('rank');
  yl(2)=ylabel('eigenvalue \lambda');
  axes(ah(4))
  xl(2)=xlabel('rank');
  axes(ah(1))
  l=legend(legs);

  longticks(ah)
  serre(ah(1:2),1/2,'across')
  serre(ah(3:4),1/2,'across')
  serre(ha(1:2),1/3,'down')
  serre(ha(3:4),1/3,'down')
  nolabels(ha(3:4),2)
  set(ah,'ylim',[-0.05 1.05])
  set([xl yl tt ah],'FontS',13)
  set([l],'FontS',12)
  set(pv,'MarkerS',4)
  movev(tt,-.02)
 case 3
  % Australia
  XY=eval(sprintf('%s(10)','australia'));
  XY(:,1)=XY(:,1)-360;
  [ah,ha,bh,th]=plotstuff('australia',XY,'ll',[],40,-1,12);
  serre(ah(1:3),4.5/3,'across')
  serre(ah(4:6),4.5/3,'across')
  serre(ah(7:9),4.5/3,'across')
  serre(ah(10:12),4.5/3,'across')
  set(ah,'xtick',[100:20:160],'xtickl',[100:20:160],...
	 'ytick',[-50:15:0],'ytickl',[-50:15:0])
  nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
 case 4
  % North America
  XY=eval(sprintf('%s(10)','namerica'));
  [ah,ha,bh,th]=plotstuff('namerica',XY,'ll',[],25.5,-1,13);
  serre(ah(1:3),1/3,'across')
  serre(ah(4:6),1/3,'across')
  serre(ah(7:9),1/3,'across')
  serre(ah(10:12),1/3,'across')
  set(ah,'xtick',[170:40:330],'xtickl',[170:40:330],...
	 'ytick',[0:30:90],'ytickl',[0:30:90])
  nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
 case 5
  % Africa
  XY=eval(sprintf('%s(10)','africa'));
  XY(:,1)=XY(:,1)-360;
  [ah,ha,bh,th]=plotstuff('africa',XY,'ll','yes',[],[],[],fi);
  if fi==0
    serre(ah(1:3),10/3,'across')
    serre(ah(4:6),10/3,'across')
    serre(ah(7:9),10/3,'across')
    serre(ah(10:12),10/3,'across')
  end
  set(ah,'xtick',[-30:30:60],'xtickl',[-30:30:60],...
	 'ytick',[-60:30:60],'ytickl',[-60:30:60])
  nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
 case 6
  % Eurasia
  XY=eval(sprintf('%s(10)','eurasia'));
  XY(:,1)=XY(:,1)-360;
  [ah,ha,bh,th]=plotstuff('eurasia',XY,'lr','maybe',20,[],[],[],L);
  set(ah,'CameraV',6)
  serre(ha(1:4),3.25,'down')
  serre(ha(5:8),3.25,'down')
  serre(ha(9:12),3.25,'down')  
  set(ah,'xtick',[-45:45:220],'xtickl',[-45:45:220],...
	 'ytick',[0:30:90],'ytickl',[0:30:90])
  nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
 case 7
  % South America
  XY=eval(sprintf('%s(10)','samerica'));
  [ah,ha,bh,th]=plotstuff('samerica',XY,'ll');
  serre(ah(1:3),12.5/3,'across')
  serre(ah(4:6),12.5/3,'across')
  serre(ah(7:9),12.5/3,'across')
  serre(ah(10:12),12.5/3,'across')
  set(ah,'xtick',[270:30:330],'xtickl',[270:30:330],...
	 'ytick',[-80:25:30],'ytickl',[-80:25:30])
  nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
 case 8
  % Amazon basin
  XY=eval(sprintf('%s(10)','amazon'));
  fi=1;
  [ah,ha,bh,th]=plotstuff('amazon',XY,'ll',[],[],[],[],fi,25);
  set(ah,'xtick',[270:20:320],'xtickl',[270:20:320],...
	 'ytick',[-30:10:10],'ytickl',[-30:10:10])

  if fi==0
    serre(ah(1:3),8/3,'across')
    serre(ah(4:6),8/3,'across')
    serre(ah(7:9),8/3,'across')
    serre(ah(10:12),8/3,'across')
    nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
  else
    serre(ah(1:4),1/2,'across')
    serre(ah(5:8),1/2,'across')
    serre(ah(9:12),1/2,'across')
    serre(ha(1:3),1,'down')
    serre(ha(4:6),1,'down')
    serre(ha(7:9),1,'down')
    serre(ha(10:12),1,'down')
    nolabels(ah(1:8),1); nolabels(ha(4:12),2); deggies(ah)
  end
 case 9
  % Orinoco basin
  XY=eval(sprintf('%s(10)','orinoco'));
  [ah,ha,bh,th]=plotstuff('orinoco',XY,'ll',[],[],[],[],[],25);
  serre(ah(1:3),2/3,'across')
  serre(ah(4:6),2/3,'across')
  serre(ah(7:9),2/3,'across')
  serre(ah(10:12),2/3,'across')
  set(ah,'xtick',[270:20:320],'xtickl',[270:20:320],...
	 'ytick',[-30:10:10],'ytickl',[-30:10:10])
  nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
 case 10
  % Greenland
  XY=eval(sprintf('%s(10)','greenland'));
  [ah,ha,bh,th]=plotstuff('greenland',XY,'ll',[],[],[],9,1,45);
  set(ah,'xlim',[262 360])
  set(ah,'ylim',[50 90])
  set(ah,'xtick',[270:30:360],'xtickl',[270:30:360],...
	 'ytick',[50:10:90],'ytickl',[50:10:90])
  nolabels(ah(1:8),1); nolabels(ha(4:12),2); deggies(ah)
  for ind=1:3
    serre(ah([1:4]+(ind-1)*4),1/3,'across')
  end
  movev(ah(1:4),-0.175)  
  movev(ah(9:12),0.175)  
  set(ah,'camerav',9.75)
 case 11
  % GPS North America
  [ah,ha,bh,th,H]=plotstuff('gpsnamerica',[],'ll',[],[],[],9,1,L);
  set(ah,'xlim',[230 252])
  set(ah,'ylim',[27 48])
  set(ah,'xtick',[230:5:250],'xtickl',[230:5:250],...
	 'ytick',[30:5:45],'ytickl',[30:5:45])
  nolabels(ah(1:8),1); nolabels(ha(4:12),2); deggies(ah)
  serre(H,1/3,'across')
  serre(H',1/3,'down')
  set(ah,'camerav',8)  
 case 12
  % Antarctica - keep editing below
  [ah,ha,bh,th,H]=plotstuff('antarctica',[],'ll',[],15,[],9,fi,L);
  % Track the amount of counterrotation would be needed
  [XY,lonc,latc]=antarctica;
  set(ah,'xlim',[15 80])
  set(ah,'ylim',[-30 20])
  % Note that I am adjusting the ticks, not doing the rotation
  set(ah,'xtick',[15:15:80],'xtickl',[15:15:80]+round(lonc),...
	 'ytick',[-30:10:20],'ytickl',[-30:10:20]+ceil(latc))
  if fi==1
    nolabels(ah(1:8),1); nolabels(ha(4:12),2); deggies(ah)
    serre(H,1/3,'across')
    serre(H',1/3,'down')
    set(ah,'camerav',8.5)
    movev(ah(1:4),-.075*2)
    movev(ah(5:8),-.075)
  else
    nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
    serre(H,1/3,'across')
    movev(ah(1:3),.01*3)
    movev(ah(4:6),.01*2)
    movev(ah(7:9),.01)
    moveh(ha(1:4),.05*2)
    moveh(ha(5:8),.05)
  end  
 case 13
  % The world's oceans
  [ah,ha,bh,th,H]=plotstuff('alloceans',[],'ul',[],15,[],9,fi,L);
  set(ah,'xtick',[0:60:360],'xtickl',[0:60:360],...
	 'ytick',[-90:45:90],'ytickl',[-90:45:90])
  nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
  axes(ah(11))
  % Fake the whole colorbar, see inside PLOTSTUFF
  colormap('kelicol')
  cb=colorbarf('hor',8,'Helvetica',[0.4445 0.06 0.2134 0.02]);
  set(cb,'xtick',linspace(indeks(get(cb,'xlim'),1),...
			  indeks(get(cb,'xlim'),2),3),...
	 'xtickl',{'-max(abs)','0','max(abs)'})
  longticks(cb,2)
 otherwise
 error('Specify valid region')
end

% Prepare outputs
varns={ah,ha,bh,th};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah,ha,bh,th,H]=plotstuff(reg,XY,lox,swti,num,swl,fozo,fi,L)
% Number of basis function to show
defval('np',12)
% Legend location
defval('lox','ul')
% Tick mark switching option
defval('swti','no')
% Axis expansion
defval('num',40)
% Data exaggeration
defval('swl',1)
% Font Size
defval('fozo',11)
% Panel geometry
defval('fi',0);
% Bandwidth
defval('L',18)

if fi==0
  [ah,ha,H]=krijetem(subnum(4,3));
else
  [ah,ha,H]=krijetem(subnum(3,4));
end
[dems,dels]=addmon(L);

% Load the required number of eigenfunctions - more for alloceans 
infx=30;
infl=[1+(infx-1)*strcmp(reg,'alloceans')];
[V,C,jk1,jk2,XYZ]=localization(L,reg,[],np*infl);
defval('XY',XYZ)

% If it's all oceans, plot every  fifth eigenfunction or so

% Modify to do only partial reconstruction to save time
r=repmat(NaN,[181 361 np]);
if L>48
  h=waitbar(0,sprintf(...
      'Expanding spherical harmonics up to degree %i',L));
end
% Get a file with the harmonics to reuse
[r(:,:,1),lor,lar,Plm]=plm2xyz([dels dems C{1}],1);
for index=1:np
  whichone=1+(index-1)*infl;
  if L>48
    waitbar(index/12,h)
  end
  if index>1
    % In the next line should expand only in the region to save time
    r(:,:,index)=plm2xyz([dels dems C{whichone}],1,[],[],[],Plm);
  end
  ka=swl*r(:,:,index);
  % Use setnans for this, rather
  ka=ka/max(max(abs(ka)));
  ka(abs(ka)<0.01)=NaN;
  axes(ah(index))
  sax=[-1 1];
  cola='kelicol';
  if strcmp(swti,'no')
    imagefnan([0 90-100*eps],[360 -90],ka,cola,sax)
    xtix=[0:90:360];
  elseif strcmp(swti,'yes')
    ka=ka(:,[181:361 1:180]);
    imagefnan([-180 90],[180 -90],ka,cola,sax);
    xtix=[-180:90:180];
  elseif strcmp(swti,'maybe')
    ka=ka(:,[231:361 1:230]);
    imagefnan([-130 90],[230 -90],ka,cola,sax);
    xtix=[-130:90:230];
  end
  ytix=[-90:45:90];

  axis image
  set(ah,'FontS',fozo-2)

  if ~strcmp(reg,'alloceans')
    hold on
    pc=plot(XY(:,1),XY(:,2),'k');
    axis(xpand([minmax(XY(:,1)) minmax(XY(:,2))],num))
    hold off
    t{index}=sprintf('%s =%7.3f','\lambda',V(index));
  else
    [axl,pc]=plotcont;
    axis([0 360 -90 90])
    set(pc,'Linew',0.5);
    %t{index}=sprintf('%s_{%i} =%9.6e','\lambda',whichone,V(whichone));
    t{index}=sprintf('%s = %i','\alpha',whichone);
    title(sprintf('%s = %.13g','\lambda',V(whichone)));
  end
  % Box labeling
  [bh(index),th(index)]=boxtex(lox,ah(index),t{index},fozo,[],[],1.1,0.8,1.2);
end
if L>48
  close(h)
end
set(th,'FontS',fozo-1)
longticks(ah)
set(ah,'xgrid','off','ygrid','off')

if fi==0
  nolabels(ah(1:9),1)
  nolabels(ha(5:12),2)
  serre(ah(1:3),1/2,'across')
  serre(ah(4:6),1/2,'across')
  serre(ah(7:9),1/2,'across')
  serre(ah(10:12),1/2,'across')
  serre(ha(1:4),2/3,'down')
  serre(ha(5:8),2/3,'down')
  serre(ha(9:12),2/3,'down')
end
seemax(ah,3) 
fig2print(gcf,'landscape')  
figdisp('sdwregions',sprintf('%s_%i',reg,L))
