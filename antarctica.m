function varargout=antarctica(res,buf,rotb)
% [XY,lonc,latc]=antarctica(res,buf,rotb)
% ANTARCTICA(...) % Only makes a plot
%
% Finds the coordinates of Antarctica, but rotates them to an equatorial
% location so KERNELC doesn't choke on the calculation.
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% buf      Distance in degrees that the region outline will be enlarged
%          by BUFFERM, not necessarily integer, possibly negative
%          [default: 0]
% rotb     0 You want the coordinates on the equator for mathematical 
%            operations (e.g. you need to make a kernel) [default]
%          1 You want the coordinates rotated back to their original
%            location (e.g. for plotting)
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
% lonc     The amount by which you need to rotate it back over z
% latc     The amount by which you need to rotate it back over y
%
% See also PLM2ROT, GEOBOXCAP, KLMLMP2ROT
% 
% Last modified by charig@princeton.edu, 03/29/2012

defval('res',0)
defval('buf',0)
defval('rotb',0)


% First part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c11=[0 -62];
cmn=[360 -83];
    
% Do it! Make it, load it, save it
XY=regselect('antarctica',c11,cmn,[],res,buf);
     
% Do some extra here to get the lonc latc
    
[~,~,XY2]=plotcont(c11(:),cmn(:),[],0);
close
% Get rid of common NaNs
XY2=XY2(~isnan(XY2(:,1)) & ~isnan(XY2(:,2)),:);
% Get rid of common NaNs
XY2=XY2(~isnan(XY2(:,1)) & ~isnan(XY2(:,2)),:);
% Find the geographical center and the area
[lonc,latc,A]=rcenter([XY2(:,1) XY2(:,2)]);
lonc=-45;


if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
  axis([-30-lonc 30-lonc -30 20])
end

% Do we return rotated coordinates?
if rotb==1
   [thetap,phip,rotmats]=rottp((90-XY(:,2))*pi/180,XY(:,1)/180*pi,-lonc*pi/180,latc*pi/180,0);
   XY = [phip*180/pi 90-thetap*180/pi];
end

% Make sure the coordinates make sense
XY(:,1)=XY(:,1)-360*[XY(:,1)>360];
XY(:,1)=XY(:,1)+360*[XY(:,1)<0];
        
% Prepare Output
varns={XY,lonc,latc};
varargout=varns(1:nargout);



