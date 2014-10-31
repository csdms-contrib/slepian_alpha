function varargout=plotprops(x,y,z,bnd,axl,ms,mrk,clox,cmap)
% PLOTPROPS(x,y,z,bnd,axl,ms,mrk,clox,cmap)
% [p,cb]=PLOTPROPS(...)
%
% Plots the third dimension as color-coded symbols on a map.
%
% INPUT:
%
% x            x-coordinate
% y            y-coordinate
% dep          z-coordinate
% bnd          Color boundaries [default: 11 for 10 bins; movable edges]
%              OR: a single number (>=4) of bin boundaries [default: 11]
% axl          Axis limits
% ms           Marker size [default: 6]
% mrk          Marker symbol [default: 'o']
% clox         Colorbar location [default: 'hor']
% cmap         Colormap string [default: grey scale]
%
% OUTPUT:
%
% p            Plot symbol handles
% cb           Color bar handle
%
% EXAMPLE:
%
% plotprops(linspace(100,200,12),repmat(0,1,12),...
%          [1:12]+0.1,1:10,[90 210 -10 10],15,[],[],'kelicol')
%
% SEE ALSO: PLOTEQ, PLOTPROP
%
% Last modified by fjsimons-at-alum.mit.edu, 06/28/2010

% Supply defaults
defval('ms',6)
defval('mrk','o')
defval('axl',[min(x)-range(x)/10 max(x)+range(x)/10 ...
	      min(y)-range(y)/10 max(y)+range(y)/10]);
defval('bnd',linspace(min(z),max(z),11));
if prod(size(bnd))==1
   bnd=linspace(min(z),max(z),bnd);
end
defval('clox','hor');
defval('flag',0);

% Plot all the data as separate symbols
for index=1:length(x)
  p(index)=plot3(x(index),y(index),index,mrk);
  hold on
end
view(2)
axis(axl)

% The following is adjusted from earlier versions, PLOTPROP and PLOTEQ
defval('cmap',flipud(gray(length(bnd)-1)))
if isstr(cmap)
  % Note that this does NOT work for KELICOL
  cmap=eval(sprintf('%s(%i)',cmap,length(bnd)-1));
  flag=1;
end

% There are length(bnd) bin starting points
zco=repmat(NaN,size(z));

% Use inequalities & find out in which bin the data are
% The first and the last bin contain the under- and overflow and so
% the first and last bin boundaries are movable in a sense; the color
% indices now run from length(bnd) down to 2, i.e. length(bnd)-1 bins
zco(z<bnd(2))=length(bnd);
for index=2:length(bnd)-1
  zco(z>=bnd(index))=length(bnd)-index+1;
end

if flag==1
  % This is the nearest-neighbor index into the color map
  zco=round(scale(zco,[1 length(bnd)-1]));
else
  % The grey scale, of course, is from 0 to 1
  zco=scale(zco,[0 10]);
end
  
% Assign the proper colors to the symbols
for index=1:length(x)
  if flag==1
    zcol=cmap(zco(index),:);
  else
    zcol=grey(zco(index));
  end
  set(p(index),'MarkerFaceColor',zcol,...
	       'MarkerEdgeColor','k',...
	       'MarkerSize',ms)
end
box on

% Add the color bar and mark all but the edges
cb=colorbar(clox);
% This is the UNIQUE colormap for the whole figure, and this it
% doesn't work with multiple panels having different color maps
colormap(cmap)

switch clox
 case 'hor'
    lims='xlim'; tix='Xtick'; tixl='Xtickl';
 case 'ver'
    lims='ylim'; tix='Ytick'; tixl='Ytickl';	
 otherwise
  error('Specify a legal colorbar orientation')
end

% Set the boundaries without further changes; feed the program
% reasonable ones... Here are the divisions
% divs=linspace(indeks(get(cb,lims),1),...
% indeks(get(cb,lims),2),length(bnd));
divs=linspace(1,length(bnd),length(bnd));
% If all data are contained within the bin boundaries, also show the
% end points; not if they don't
if min(z)>=bnd(1);   a=1; else a=2; end
if max(z)<=bnd(end); b=0; else b=1; end

% Here you piece them together from the limits - this is best
set(cb,tix,divs(a:end-b),tixl,bnd(a:end-b))
% Here you get them from whatever Matlab had already decided
% gotem=get(cb,tix);
% set(cb,tix,gotem(a:end-b),tixl,bnd(a:end-b))
longticks(cb,2)

varns={p,cb};
varargout=varns(1:nargout);

hold off

