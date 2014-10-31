function [imagemat,minmat,maxmat]=dat2col(matrix,seeaxis,cmapnum,varargin)
% [imagemat,minmat,maxmat]=DAT2COL(matrix,seeaxis,cmapnum,cmlenth)
%
% Maps the 'matrix' onto a colormap of which you give the
% number 'cmapnum' in a sequence of multiple colormaps
% that are formed as CMAP=[cmap1 ; cmap2 ; ... ; cmap3]
% where the individual cmapx are 64*3, specifying also the 
% 'seeaxis' that define the extreme values of the colormap.
% This then either means data truncation to fit within 'seeaxis'
% or stretching - just as 'caxis' would work. After this,
% all you have to do is 'image', not 'imagesc' since the explicit
% indices in the master color map are returned.
% Also returns the original min and max of the data, in case you need
% them for axis labeling. Default 'seeaxis' is all the data.
% Default 'cmlenth' is 64.
%
% EXAMPLE:
%
% clear all ; load clown
% cmap1=jet; cmap2=gray(size(jet,1)); cmap3=copper; cmap4=autumn;
% cmap=[cmap1 ; cmap2 ; cmap3 ; cmap4];
% im1=dat2col(X,[1 80],1); im2=dat2col(X,[1 60],2); 
% im3=dat2col(X,[1 60],3); im4=dat2col(X,[-60 60],4); 
% subplot(221) ; image(im1) ; subplot(222) ; image(im2)
% subplot(223) ; image(im3) ; subplot(224) ; image(im4)
% colormap(cmap)
%
% Last modified by fjsimons-at-alum.mit.edu, 06/07/2007

% Default colormap length
if nargin~=4
  deflen=64;
else
  deflen=varargin{1};
end
%------------------------
[m,n]=size(matrix);
minmat=min(matrix(:));
maxmat=max(matrix(:));

if isempty(seeaxis)
  seeaxis=[minmat maxmat];
end
 
seeaxis=sort(seeaxis);

% Figure out what the colormap index limits are  going to be
cmlo=1+(cmapnum-1)*deflen;
cmhi=deflen+(cmapnum-1)*deflen;

% Saturate - if need be
matrix(matrix<=seeaxis(1))=seeaxis(1);
matrix(matrix>=seeaxis(2))=seeaxis(2);

% Strecth - if need be
matrix=[matrix(:) ; seeaxis(1) ; seeaxis(2)];

% Scale to appropriate colormap and put back in shape
matrix=scale(matrix,[cmlo cmhi]);
imagemat=reshape(matrix(1:end-2),m,n);
