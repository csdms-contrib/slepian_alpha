function varargout=plottopo(data,contcut,caxcon,caxoc,c11,cmn)
% H=PLOTTOPO(data,contcut,caxcon,caxoc,c11,cmn)
%
% Plots "topography" on an image map, using good colormap.
%
% INPUT:
%
% data           The data matrix, where (1,1) is NorthWest
% contcut        Switch color maps at this elevation
% caxcon         The caxis limits for the continental color map
% caxoc          The caxis limits for the oceanic color map
% c11            The coordinates of the pixel center for (1,1)
% cmn            The coordinates of the pixel center for (M,N)
%
% OUTPUT:
%
% H              An axis handle to the image map
%
% See also PLOTGRAV, ADDCB, CAX2DEM, JOINCOLMAP
%
% Last modified by fjsimons-at-alum.mit.edu, 10/08/2008

defval('c11',[0 90])
defval('cmn',[360 -90])

[dem,dax,ziro]=sergeicol;
defval('contcut',0);
defval('caxcon',[0 1500]);
defval('caxoc',[-7000 0]);

Tcont=find(data>=contcut);
Toc=find(data<contcut);

% Make the actual colormap from the two portions of sergeicol
RGB=joincolmap(data,[],Tcont,dem(ziro+1:end,:),caxcon);
RGB=joincolmap(data,RGB,Toc,dem(1:ziro,:),caxoc);

% Index this color map directly
if isempty(c11)
  h=image(RGB);
else
  h=imagefdir(c11,cmn,RGB);  
end

warning off
axis image
warning on

if nargout
  varargout{1}=h;
end

