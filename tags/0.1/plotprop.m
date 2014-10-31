function varargout=plotprop(lon,lat,prop,bnd,axl,ms,mrk,clox)
% PLOTPROP(lon,lat,prop,bnd,axl,ms,mrk,clox)
% [ah,cb,p,pc,pp]=PLOTPROP(...)
%
% Plots the third dimension as color-coded symbols on a map.
%
% INPUT:
%
% lon          Longitude (x-coordinate)
% lat          Latitude (y-coordinate)
% dep          Depth or some other property
% bnd          Color boundaries [default: 10 bins]
% axl          Axis limits
% ms           Marker size [default: 6]
% mrk          Marker symbol [default: 'o']
% clox         Colorbar location [default: 'hor']
%
% OUTPUT:
%
% ah           Figure panel handles           
% cb           Color bar handle
% p            Plot symbol handles
% pc           Continental outline handles
% pp           Plate boundary handles
%
% EXAMPLE:
%
% plotprop(linspace(100,200,12),repmat(0,1,12),...
%          [1:12]+0.1,1:10,[90 210 -10 10],15)
% plotprop(linspace(100,200,12),repmat(0,1,12),...
%          [0:11]+0.1,1:10,[90 210 -10 10],15)
% plotprop(linspace(100,200,12),repmat(0,1,12),...
%          [-2:9]+0.1,1:10,[90 210 -10 10],15)
% plotprop(linspace(100,200,12),repmat(0,1,12),...
%          [1:12]+0.1,1:2:20,[90 210 -10 10],15)
%
% SEE ALSO: PLOTEQ
%
% Last modified by fjsimons-at-alum.mit.edu, 07/07/2008

% Supply defaults
defval('ms',6)
defval('mrk','o')
defval('axl',[-120+360 -115+360 32 37]);

defval('bnd',linspace(min(prop),max(prop),10));

defval('clox','hor');

lon=lon+360*(lon<0);

% California
ah=gca;
for index=1:length(lon)
  p(index)=plot3(lon(index),lat(index),index,mrk);
  hold on
end

view(2)
[a,b]=plotcont;
c=plotplates;

axis image
axis(axl)

% There are length(bnd) bin starting points
propco=repmat(NaN,size(prop));
if min(prop)<=bnd(1)
  %error('Data minimum lower than bin minimum - adjust bins')
  bnd(1)=floor(min(prop))-eps;
end
% Use inequalities & find out in which bin the data are
for index=1:length(bnd)
  propco(prop>bnd(index))=length(bnd)-index+1;
end

% The grey scale, of course, is from 0 to 1
% always include the lowest, too
propco=indeks(scale([propco(:)' 1],[0 10]),'1:end-1');

% Assign the proper colors to the symbols
for index=1:length(lon)
  set(p(index),'MarkerFaceColor',grey(propco(index)),...
	       'MarkerEdgeColor','k',...
	       'MarkerSize',ms)
end
box on

cb=colorbar(clox);
colormap(flipud(gray(length(bnd))))
caxis([0 1])

for index=1:length(bnd)
  bon{index}=sprintf('>%i',round(bnd(index)));
end

switch clox
 case 'hor'
  set(cb,'XTick',1.5:length(bnd)+0.5,'XTickL',bon)
 case 'ver'
  set(cb,'YTick',1.5:length(bnd)+0.5,'YTickL',bon)
  otherwise
  error('Specify a legal colorbar orientation')
end

set(cb,'TickDir','out','TickLength',[0.02 0.025]/3)

longticks(ah,2)

varns={ah,cb,p,b,c};
varargout=varns(1:nargout);

hold off



