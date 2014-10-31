function varargout=addcb(pos,caxcon,caxoc,parm,tint,invt)
% [cb,xcb]=addcb(pos,caxcon,caxoc,parm,tint,invt)
%
% Adds an explicit colorbar - suitable for use with multiple or
% complicated color maps, or with directly indexed RGB maps, as in
% IMAGEFNAN. Remember, only the tick LABELS are meaningful; their actual
% positions are not. 
%
% INPUT:
%
% pos         'hor' [default], 'vert' or [lbwh]
% caxcon      Caxis limits for the "continental" colormap
% caxoc       Caxis limits for the "oceanic" colormap
%             The defaults of the above are for a dual colormap
%             suitable for global topography; if another colormap is
%             specified, min([caxoc caxcon]) and max([caxoc caxcon])
%             are used as the axis limits on the colorbar; so just give
%             them the same values if you're using only one as with IMAGEFNAN.
% parm        1 Topography colormap [default]
%             2 Gravity colormap composed of two sections
%             3 Another colormap identified by a string, e.g. 'kelicol'
%             4 An actual colormap as a rgb matrix, e.g. gray(21)
% tint        Interval between the tickmarks [default: one tenth of the range]
% invt        1 Invert the color bar in question
%             0 Don't [default]
%
% OUTPUT:
%
% cb           Colorbar axis handle
% xcb          Colorbar xlabel axis handle
%
% CAVEAT:
%
% When in doubt, specify the colorbar by name string
% 
% SEE ALSO:
%
% PLOTTOPO, PLOTGRAV, JOINCOLMAP, CAX2DEM, SERGEICOL, DEMMAP, IMAGEFDIR
%
% Last modified by fjsimons-at-alum.mit.edu, 10/10/2011

defval('caxcon',[0 1500]);
defval('caxoc',[-7000 0]);
defval('parm',1)
defval('pos','hor')
defval('invt',0)

ah=gca;
fpos=getpos(ah);

poso=pos;
if isstr(pos)
  cb=colorbar(pos);
  npos=getpos(ah);
  pos=getpos(cb);
  delete(cb)
  set(ah,'position',npos)
end
cb=axes('position',pos);

miC=min([caxoc caxcon]);
maC=max([caxoc caxcon]);
cbd=linspace(miC,maC,500);

defval('tint',(maC-miC)/10)

wis=2;
if length(parm)==1 && parm==1
  % Not quite right if it turns out to be vertical 
  % New length addition to skip for BW
  plottopo(cbd,[],caxcon,caxoc,[0 1],[1 0])
elseif length(parm)==1 && parm==2
  % New length addition to skip for BW
  plotgrav(cbd,-sign(cbd),0,caxcon,caxoc)
elseif [~isstr(parm) && strcmp(poso,'hor')] || ...
      [~isstr(parm) && poso(3)>poso(4)]
  % New if for BW color bars called by gray(N)
  h=parm; if invt==1; h=flipud(h); end
  h=reshape(h,[1 size(h,1) size(h,2)]);
  imagefdir([miC 1],[maC 0],h);
elseif [~isstr(parm) && strcmp(poso,'vert')]  || ...
      [~isstr(parm) && poso(4)>poso(3)]
  % New if for BW color bars called by gray(N)
  h=parm; if invt==0; h=flipud(h); end
  h=reshape(h,[1 size(h,1) size(h,2)]);
  for index=1:3
    hh(:,:,index)=h(:,:,index)';
  end
  imagefdir([miC 1],[maC 0],hh);
  wis=1;
elseif [isstr(parm) && strcmp(poso,'hor')] || poso(3)>poso(4)
  h=eval(parm); if invt==1; h=flipud(h); end
  h=reshape(h,[1 size(h,1) size(h,2)]);
  imagefdir([miC 1],[maC 0],h);
elseif [isstr(parm) && strcmp(poso,'vert')] || poso(4)>poso(3)
  h=eval(parm); if invt==0; h=flipud(h); end
  h=reshape(h,[1 size(h,1) size(h,2)]);
  for index=1:3
    hh(:,:,index)=h(:,:,index)';
  end
  imagefdir([miC 1],[maC 0],hh);
  wis=1;
end

noticks(cb,wis)
longticks(cb)

if [~isstr(parm) && strcmp(poso,'hor')]
  % New if for BW
  xtcb=get(cb,'xtick'); 
  xlcb=get(cb,'xtickl');
  set(cb,'xtick',xtcb(1:2:end))
  set(cb,'xtickl',cellstr(xlcb(1:2:end,:)))
elseif [~isstr(parm) && strcmp(poso,'hor')]
  % New if for BW
  ytcb=get(cb,'ytick'); 
  ylcb=get(cb,'ytickl');
  set(cb,'ytick',ytcb(1:2:end))
  set(cb,'ytickl',cellstr(ylcb(1:2:end,:)))
elseif [length(parm)==1 && parm==1] || [length(parm)==1 && parm==2] || ...
	([isstr(parm) && strcmp(poso,'hor')] || poso(3)>poso(4))
  % New length addition to skip for BW
  xtcb=get(cb,'xtick'); 
  xlcb=get(cb,'xtickl');
  set(cb,'xtick',xtcb(1:2:end))
  set(cb,'xtickl',cellstr(xlcb(1:2:end,:)))
elseif [isstr(parm) & strcmp(poso,'vert') || poso(4)>poso(3)]
  ytcb=get(cb,'ytick'); 
  ylcb=get(cb,'ytickl');
  set(cb,'ytick',ytcb(1:2:end))
  set(cb,'ytickl',cellstr(ylcb(1:2:end,:)))
end

if ~isstr(parm) && length(parm)==1
  labn={'Topography (km)' 'Gravity (mGal)'};
  xcb=xlabel(labn{parm});
elseif [~isstr(parm) && strcmp(poso,'hor')]
  xcb=xlabel('add y-label here');
elseif [~isstr(parm) && strcmp(poso,'vert')]
  xcb=ylabel('add x-label here');
elseif [isstr(parm) && strcmp(poso,'vert')] ...
      || [~isstr(poso) && poso(4)>poso(3)]
  xcb=ylabel('add y-label here');
elseif [isstr(parm) & strcmp(poso,'hor')] ...
      || [~isstr(poso) && poso(3)>poso(4)]
  xcb=xlabel('add x-label here');
else
  error('Specify valid options')
end

axes(cb)
axis normal

% Add the tick marks on the color bar
if ~isstr(poso) && poso(3)>poso(4)
  poso='hor';
elseif ~isstr(poso) && poso(4)>poso(3)
  poso='vert';
end

% Transform the tick marks
cbarticks(cb,[miC maC],tint,poso)

% Output
varns={cb,xcb};
varargout=varns(1:nargout);
