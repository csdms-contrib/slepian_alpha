function varargout=dphregion(th,N,region) 
% [phint,thp,php]=DPHREGION(th,N,region)
%
% Finds the longitude crossings at a given colatitude of the coastlines
% of some region that is to be treated as occupying flat Cartesian
% geometry (which is, admittedly, a bit of a limitation).
%
% INPUT:
%
% th         Colatitudes where you want it evaluated [degrees]
% N          Smoothness of the spline version [only for region strings]
% region     A string: 'england', 'africa', 'namerica', 'greenland',
%                       'australia', 'eurasia', etc. OR
%            A matrix with XY/[lon,lat] coordinates [degrees]
%            note that this can be several closed curves separated by
%            a row of [NaN NaN]
%
% OUTPUT:
%
% phint      Longitudinal interval(s) defining the integration domain
% thp        Colatitude matrix for plotting
% php        Longitude matrix for plotting
%
% EXAMPLE:
%
% dphregion('demo1',[],'africa')
% dphregion('demo1',[],'contshelves')
% dphregion('demo2',[],'eurasia')
% dphregion('demo2',10)
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

if ~isstr(th)
  % First, get the coordinates of region; these are sorted
  defval('N',10);
  if isstr(region)
    XY=eval(sprintf('%s(%i)',region,N));
  else
    XY=region;
  end

  % Input now in colatitude
  XY(:,2)=90-XY(:,2);
  
  % Calculate the crossings of XY at th
  [phint,thp,php]=phicurve([XY(:,2) XY(:,1)],th);
  
  % Distribute output
  varns={phint,thp,php};
  varargout=varns(1:nargout);
  % AFTER THIS NOTHING BUT DEMOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(th,'demo1')
  defval('N',10)
  defval('region','africa')
  eval(sprintf('%s(%i)',region,N));
  eval(sprintf('XY=%s(%i);',region,N));
  thN=90-max(XY(:,2));
  thS=90-min(XY(:,2));
  XY(:,2)=90-XY(:,2);
  TH=thN+rand*(thS-thN);
  hold on
  XL=xlim;
  [phint,thp,php]=phicurve([XY(:,2) XY(:,1)],TH);
  pl=plot(php,90-thp,'g-','LineW',2);
  plot(XL,[90-TH 90-TH],'k--')
  p=plot(phint,90-TH,'o');
  set(p,'MarkerE','k','MarkerF','y')
  t=title(num2str(length(phint)),'FontS',30);
  hold off
  if nargout==1
    phint=dphi;
  end  
elseif strcmp(th,'demo2')
  % Plot this for the time being
  defval('region','africa')
  defval('N',200)
  ah(2)=subplot(212);
  ah(1)=subplot(211); 
  Nk=100;
  eval(sprintf('%s(%i)',region,Nk));
  eval(sprintf('XY=%s(%i);',region,Nk));
  thN=90-max(XY(:,2));
  thS=90-min(XY(:,2));
  XY(:,2)=90-XY(:,2);
  thetas=linspace(thN,thS,N);
  hold on; p=[];
  XL=xlim;
  for ind=1:length(thetas)
    TH=thetas(ind);
    l=plot(XL,[90-TH 90-TH],'k-');
    % Add eps so it is never zero, essentially
    xmth=XY(:,2)-TH;
    dsx=diff(sign(xmth));
    % Now it can be the one before, or after the crossing, how about
    colf=find(dsx);
    % This now returns the one on the negative side of the line
    colx=colf+(dsx(colf)==-2);
    colx2=colx-(dsx(colf)==-2)+(dsx(colf)==2);
    L(ind)=length(colx);
    if mod(L(ind),2)
      warning(sprintf('Cannot find pairs of crossings at th=%8.3f',TH))
    end
    axes(ah(2))
    plot(XY(:,1),xmth,'-+'); hold on; grid on
    xls=[min(XY(colx,1))-range(XY(colx,1))/3 ...
	 max(XY(colx,1))+range(XY(colx,1))/3];
    % We want the interpolated point at which it becomes exactly zero
    xlim(xls)
    plot(xls,[0 0],'k')
    ylim([-1 1]/100)
    % Then one point was exactly hit, this is the thN case
    if all(colx==colx2)
      dph=XY([colx(2) colx2(2)]); % This works
    else
      for ond=1:L(ind)
	dph(ond)=interp1(xmth([colx(ond) colx2(ond)]),...
			 XY([colx(ond) colx2(ond)],1),0,'linear');
      end
    end
    plot(XY(colx,1),xmth(colx),'bs');
    plot(XY(colx2,1),xmth(colx2),'rv'); 
    plot(dph,repmat(0,size(dph)),'g*'); 
    clear dph
    hold off
    for indo=1:length(colx)
      axes(ah(1))
      p(indo)=plot(XY(colx(indo),1),90-XY(colx(indo),2),'o');
    end
    if ~isempty(p)
      set(p,'MarkerE','k','MarkerF','y')
      % t=text(362,58,num2str(L(ind)),'FontS',30);
      t=title(num2str(L(ind)),'FontS',30);
      pause
      delete([l p t]); p=[];
    end
  end
  hold off
  disp(sprintf('Number of uneven crosses %i',sum(mod(L,2))))
else
  error('Speficy valid demo')
end
