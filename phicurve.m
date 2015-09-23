function varargout=phicurve(thph,th)
% [phint,thp,php,forreal]=PHICURVE(thph,th)
%
% Finds the longitude crossings and (thus) integration domains of a
% closed curve parameterized in colatitude/longitude space at certain
% query points of colatitude. Note that the points must be given as
% running around the curve clockwise or anticlockwise. Note that the
% region is treated as occupying flat Cartesian geometry (which is,
% admittedly, a bit of a limitation). Could be curves separated by NaNs.
%
% INPUT:
%
% thph            Colatitude/longitude of the closed curve [degrees]
% th              Colatitude at which crossings are required [degrees]
%
% OUTPUT:
%
% phint           A matrix with crossings/intervals and zeroes, of
%                 dimensions MxN where M=length(th) and N can be anything
%                 depending on the number of oscillations of the curve
% thp             Colatitude matrix for hatched plotting, if possible
% php             Longitude matrix for hatched plotting, if possible
% forreal         Indices to the ones that are for real (could be AT zero)
%
% EXAMPLE:
%
% phicurve('demo1','africa') % A geographic region
% phicurve('demo2') % A random blob
%
% SEE ALSO:
%
% DPHREGION, DPHPATCH, DPHSQUARE, CURVECHECK
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012
% Last modified by charig-at-princeton.edu, 08/14/2015

if ~isstr(thph)
  % For every th, find the relevant phint
  xmth=repmat(thph(:,1),1,length(th))-repmat(th(:)',length(thph(:,1)),1);
  dsx=diff(sign(xmth));
  if sum((dsx(:)==0)-1)==0
    error('Specify at least one colatitude within the data range')
  end
  % Now it can be the one before, or after the crossing, how about
  [colf,colj]=find(dsx);
  colr=sub2ind(size(dsx),colf,colj);
  % This now returns the one on the negative side of the line
  colx=colf+(dsx(colr)==-2);
  colx2=colx-(dsx(colr)==-2)+(dsx(colr)==2);
  L=length(colx);
  if mod(L,2)
    error(sprintf('Cannot find pairs of crossings'))
  end
  % Then one point was exactly hit, this is the thN or thS case
  if length(colx)==2 && all(colx==colx2)
    phint=thph([colx(2) colx2(2)],2);
    thp=[th th];
    php=phint;
  else
    for ond=1:L
      % In case you have a node immediate followed by a crossing
      if colx(ond)==colx2(ond)
	phint(ond)=NaN;
      else
	phint(ond)=interp1(xmth([colx(ond) colx2(ond)],colj(ond)),...
			   thph([colx(ond) colx2(ond)],2),0,'linear');
      end
    end
    % Debate whether this is useful or not wrt to node/crossing
    %    phint=phint(~isnan(phint));
    % ACTUALLY, IF THE NAN'S ARE NOT CONSECUTIVE PAIRS GET SPECIAL CASE

    % Now rearrange back to the number of requested points
    % But there could be points with more or less than 2 crossings
    % Maximum number of times a crossing is repeated
    [a,b]=degamini(colj);
    rowj=colj;
    colj=matranges(reshape([repmat(1,length(b),1) b']',length(b)*2,1))';
    pint=repmat(NaN,length(th),max(b));    
    subsi=(colj-1)*length(th)+rowj;
    pint(subsi)=phint;
    if length(b)==length(th)
      wt=0;
      thp=reshape(gamini(th,b),2,length(phint)/2);    
    else
      wt=1;
      thp=[];
    end
    % Need to sort since contour may be given in any order
    phint=sort(pint,2);
    if wt==0
      php=reshape(phint(subsi),2,length(colj)/2);
    else
      php=[];
    end
    % Make them zero so the integral doesn't do anything
    forreal=~isnan(phint);
    phint(~forreal)=0;
  end
  varns={phint,thp,php,forreal};
  varargout=varns(1:nargout);
elseif strcmp(thph,'demo1')
  defval('N',10)
  defval('th','africa')
  region=th;
  eval(sprintf('XY=%s(%i);',region,N));
  thph=[90-XY(:,2) XY(:,1)];
  Nth=ceil(rand*300);
  th=linspace(min(thph(:,1)),max(thph(:,1)),Nth);
  [phint,thp,php]=phicurve(thph,th);
  plot(php,90-thp,'k-'); hold on
  plot(thph(:,2),90-thph(:,1),'k-');
  hold off
  axis equal; grid on
  title(sprintf('Number of crossings %i',Nth))
elseif strcmp(thph,'demo2')
  [x,y]=blob(1,1);
  thph=[y(:) x(:)];
  Nth=ceil(rand*300);
  th=linspace(min(thph(:,1)),max(thph(:,1)),Nth); 
  [phint,thp,php]=phicurve(thph,th);
  plot(php,90-thp,'k-'); hold on
  plot(thph(:,2),90-thph(:,1),'k-');
  hold off
  axis equal; grid on
  title(sprintf('Number of crossings %i',Nth))
end


