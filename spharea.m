function [Ao4p,varargout]=spharea(c11,cmn) 
% [Ao4p,XY]=SPHAREA(c11,cmn)
% [Ao4p,XY]=SPHAREA(region)
% [Ao4p,XY]=SPHAREA(XY)
% [Ao4p,XY]=SPHAREA(TH,sord)
%
% Computes FRACTIONAL surface areas on the unit sphere.
%
% INPUT:
%
% c11       [x,y]/[lon,lat]-coordinates of the upper left corner
% cmn       [x,y]/[lon,lat]-coordinates of the bottom right corner
%           All coordinates are given in degrees.
%           Both inputs may be matrices of size (Mx2); OR
% region    A string name with an approved region such as 'africa', OR
%           a cell array containing the region name and buffer such as {'greenland' 0.5}
% XY        its coordinates (such as supplied from 'africa' etc), OR
% TH        A colatitudinal radius [in degrees]
% sord      0 Allowable option for region string
%           1 Single cap of radius TH [no default - must specify]
%           2 Double cap left by SUBTRACTING belt of width 2TH
%           3 Equatorial belt of width 2TH
%
% OUTPUT:
%
% Ao4p      Area as a fraction of the area of the unit sphere.
%           Multiply by 4*pi*radius^2 for a known sphere.
% XY        The outlines of the region used, if available
%
% EXAMPLE I: the fractional surface area of the Earth is
%
%      spharea([0 90],[360 -90]) % and the continental coverage
%      spharea('africa')+spharea('eurasia')+spharea('namerica')+...
%         spharea('samerica')+spharea('greenland')+spharea('australia')
% 
% EXAMPLE 2: area between one degree of longitude for all latitudes -
% this goes as the cosine of the latitude of course
% 
%      latsup=90:-1:-89; latsdwn=89:-1:-90;
%      lonslft=zeros(size(latsup)); lonsrgt=ones(size(latsup));
%      c11=[lonslft' latsup']; cmn=[lonsrgt' latsdwn']; subplot(121)
%      imagesc([0.5 0.5],[90 -90],spharea(c11,cmn)); colormap gray; axis xy
%      subplot(122); 
%      imagesc([0.5 0.5],[90 -90],cos(linspace(-pi/2,pi/2,length(latsup)))');
%      axis xy
%
% EXAMPLE 3:
%
% [a,xy]=spharea('australia'); a-areaint(xy(:,2),xy(:,1))
%
% SEE ALSO:
%
% RCENTER, AREAINT
%
% Last modified by charig-at-princeton.edu, 05/14/2015
% Last modified by fjsimons-at-alum.mit.edu, 07/21/2014

defval('XY',NaN)

if nargin==2 && ~all(cmn(:)==0)
  if all(size(c11)==size(cmn)) && length(cmn)~=1
    % Conversion to radians
    c11=c11*pi/180;
    cmn=cmn*pi/180;
    
    lon1=c11(:,1);
    lat1=c11(:,2);
    lon2=cmn(:,1);
    lat2=cmn(:,2);
  
    Ao4p=abs(lon1-lon2).*abs(sin(lat1)-sin(lat2))/4/pi;
  elseif length(cmn)==1 
    % This means that the variable SORD is specified
    TH=c11;
    sord=cmn;
    switch sord
     case 1
      % Area of the single cap of width TH
      Ao4p=(1-cos(TH*pi/180))/2;
     case 2
      % Area of the double cap each of width 90-TH
      Ao4p=(1-sin(TH*pi/180));
     case 3
      % Area of the belt of width 2TH
      Ao4p=sin(TH*pi/180);
    end    
  end
else
  defval('ngl',200)
  if isstr(c11)
    % Specify the region by name
    defval('N',10);
    region=c11;
    eval(sprintf('XY=%s(%i);',region,N));
  elseif iscell(c11)
    % Specify the region by name
    defval('N',10);
    region=c11{1};
    eval(sprintf('XY=%s(%i,%f);',region,N,c11{2}));
    % Do this so dphregion below is correct for these coordinates
    region=XY; N=NaN;
  else
    % Specify the region by the coordinates of the bounding curve
    XY=c11;
    region=XY;
    N=NaN;
  end

  % Calculate northernmost and southernmost colatitude
  thN=90-max(XY(:,2)); thN=thN*pi/180;
  thS=90-min(XY(:,2)); thS=thS*pi/180;
  % Calculate Gauss-Legendre points
  intv=cos([thS thN]);
  % Make this number high
  nGL=min(ngl,size(XY,1)/2);
  [w,x,Ngg]=gausslegendrecof(nGL,[],intv);
  % Now get the longitudinal intersections
  [phint,thp,php]=dphregion(acos(x)*180/pi,N,region);
  % plot(thp,php)
  phint=phint*pi/180;
  % This is not ridiculous as I may be giving it multiple intervals
  I=coscos(acos(x),0,0,phint);
  Ao4p=w(:)'*I/4/pi;
end

% Provide desired output
varns={XY};
varargout=varns(1:nargout-1);
