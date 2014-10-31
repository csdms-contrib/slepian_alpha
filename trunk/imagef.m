function varargout=imagef(c11,cmn,matrix,dlat)
% h=IMAGEF(c11,cmn,matrix,dlat)
%
% Uses IMAGESC to plot data in a proper geographic reference
%
% INPUT:
%
% c11         Physical coordinates of matrix element (1,1) [default:   0  90]
% cmn         Physical coordinates of matrix element (M,N) [default: 360 -90]
% matrix      The data to be plotted
% dlat        The interval for the latitude tick marks [default: 45]
%
% OUTPUT:
%
% h           A graphics handle to the plot
%
% NOTE:
%
% This is for PIXEL-CENTERED COORDINATES. Should C11 and CMN be instead
% in GRID/NODE CENTERED COORDINATES, to give the boundaries of the grid, use
% IMAGEF(c11+[dlo -dla]/2,cmn+[-dlo dla]/2,both) with DLO and DLA the X
% and Y sampling intervals. See, e.g. the ETOPO1 webpages
%
% See also: IMAGEFDIR, ADDCB, IMAGEFNAN
%
% Last modified by fjsimons-at-alum.mit.edu, 07/24/2010

if nargin==1
  h=imagef([],[],c11);
  return
end

defval('c11',[0 90])
defval('cmn',[360 -90])

cm1=[c11(1) cmn(2)];
c1n=[cmn(1) c11(2)];
h=imagesc([cm1(1) c1n(1)],[cm1(2) c1n(2)],flipud(matrix));  axis xy

% Equivalent pcolor
% pcolor(linspace(c11(1),cmn(1),size(both,2)+1),...
%    linspace(c11(2),cmn(2),size(both,1)+1),adrc(matrix));

defval('dlat',45)
if all(c11==[0 90]) & all(cmn==[360 -90])
  set(gca,'ytick',[-90:dlat:90])
  set(gca,'xtick',[0:90:360])
  deggies(gca)
end

varns={h};
varargout=varns(1:nargout);

