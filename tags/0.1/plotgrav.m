function varargout=plotgrav(grav,topo,contcut,caxcon,caxoc,c11,cmn)
% H=PLOTGRAV(grav,topo,contcut,caxcon,caxoc,c11,cmn)
%
% Plots "gravity" data in different color scheme
% according to the topography.
%
% See also PLOTTOPO
%
% Last modified by fjsimons-at-alum.mit.edu, 10/08/2008

defval('contcut',-1500);
defval('caxcon',[-350 0]);
defval('caxoc',[0 300]);
defval('c11',[])
defval('cmn',[])

% Use this colormap
colormap default ; def=colormap;

Bcont=find(topo>=contcut);
Boc=find(topo<contcut);
RGB=joincolmap(grav,[],Bcont,def,caxcon);
RGB=joincolmap(grav,RGB,Boc,flipud(gray(size(def,1))),caxoc);

if isempty(c11)
  h=image(RGB);
else
  h=imagefdir(c11,cmn,RGB);  
end

if nargout
  varargout{1}=h;
end
