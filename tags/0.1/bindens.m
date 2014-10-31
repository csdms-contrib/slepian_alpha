function [h,c11,cmn,hh,ybine]=bindens(x,y,nxbins,nybins)
% [h,c11,cmn,hh,ybine]=BINDENS(x,y,nxbins,nybins)
%
% Density plot for data in 2D histograms
%
% INPUT:
%
% x,y              The data vectors
% nxbins, nybins   The number of bins in the x and y direction
%
% OUTPUT:
%
% h          The '2D' histogram
% c11,cmn    The centers of the top left and bottom right of this histogram
% hh         Globally normalized histogram
% ybine      The y bin edges that are being used
%
% See also ROW2STATS, HIST2D
%
% Last modified by fjsimons-at-alum.mit.edu, 03/18/2013

defval('nxbins',10)
defval('nybins',10)
defval('xbin',range(x)/nxbins);
defval('ybin',range(y)/nybins);

% Specify the y-bin edges, that's one more than there are bins
ybine=linspace(min(y),max(y),nybins+1);

% Sort the data in the first column
[x,I]=sort(x,1);

% And have the second column follow
y=y(I);

% Bin the x-data
ix=ceil((x-min(x))/xbin);
ix=ix+(ix==0);

% Also must put in nans for the bins that didn't happen
adix=skip(1:max(ix),unique(ix));
ix=[ix ; adix(:)];
y=[y ; nan(size(adix(:)))]; 

% Now use ROW2STATS to get the y-histograms
[g,s,h,hh]=row2stats(ix,y,ybine);

% Do not do flipud as the histogram treats the bins as one-sided to the
% right but rather find the pixel-centered coordinates of the histogram?
c11(1)=min(x)+xbin/2;
cmn(1)=max(x)-xbin/2;
c11(2)=min(y)+ybin/2;
cmn(2)=max(y)-ybin/2;
