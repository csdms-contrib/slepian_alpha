function varargout=imagefnan(c11,cmn,matrix,colmap,cax,nancol,invo,setnan)
% [h,cax,Xrgb,c11,cmn]=IMAGEFNAN(c11,cmn,matrix,colmap,cax,nancol,invo,setnan)
% [h,cax,Xrgb,c11,cmn]=IMAGEFNAN(matrix)
% 
% Uses IMAGE to plot data in a proper geographic reference
%
% INPUT:
%
% c11         Physical coordinates of matrix element (1,1) [default: 0 1]
% cmn         Physical coordinates of matrix element (M,N) [default: 1 0]
% matrix      The data to be plotted
% colmap      A (string with the name of the) colormap [default: kelicol]
% cax         Color saturation levels [default: HALVERANGE]
% nancol      The color assigned to NaN's [default: white]
% invo        1 Invert the current colormap
%             0 Don't invert the current colormap [default]
% setnan      1 Set values smaller than 100*eps to NaN [default]
%             0 Don't 
%             A scalar that defines the threshold compared to max(abs) (a denominator!)
%
% OUTPUT:
%
% h           The axis handle to the image being plotted
% cax         The color limits being used
% Xrgb        The RGB matrix being plotted
% c11,cmn     The physical coordinates
%
% See IMAGEF, IMAGEFDIR, JOINCOLMAP, ADDCB, HALVERANGE
%
% Last modified by fjsimons-at-alum.mit.edu, 06/24/2014

if nargin==1
  % Then the variable "c11" really is the variable "matrix"
  h=imagefnan([],[],c11,[],[],[],[],1);
  return
end

defval('c11',[0 1])
defval('cmn',[1 0])
defval('colmap',kelicol)
defval('nancol',[1 1 1])
% defval('cax',round(halverange(matrix,[],NaN)))
defval('cax',halverange(matrix,[],NaN))
defval('invo',0)
defval('setnan',1)

if prod(size(matrix))>1e8
  error('You likely going into swap might not like this')
end

if setnan==1
  matrix(abs(matrix)<eps*100)=NaN;
elseif setnan~=0
  [matrix,thre]=setnans(matrix,setnan);
  disp(sprintf('Using threshold %g',thre))
end

colmaps=colmap;

if isstr(colmap)
   colmap=eval(colmap);
end
if invo==1
  colmap=flipud(colmap);
end

ncol=linspace(cax(1),cax(2),size(colmap,1))';
thenans=find(isnan(matrix));

matrix(thenans)=0;
matrix(matrix<cax(1))=cax(1);
matrix(matrix>cax(2))=cax(2);

% What if there is only one non-NaN?
if ~~sum(diff(ncol))
  % Interpolate the color map
  X=interp1(ncol,colmap,matrix(:),'nearest');
else
  % Take a single color, whatever
  X=repmat(1/pi,length(matrix(:)),3);
end

if ~isempty(thenans),
  X(thenans(:),1)=nancol(1);
  X(thenans(:),2)=nancol(2);
  X(thenans(:),3)=nancol(3);
end

X=reshape(X,[size(matrix),3]);

% This is the call in case you want to substitute specific values
% Of course you could get X from get(gcf,'CData')
%h=image([c11(1) cmn(1)],[cmn(2) c11(2)],flipdim(X,1));
h=image([c11(1) cmn(1)],[c11(2) cmn(2)],X);
axis xy image

if nargin>4
  % Suggests how to add a simple color bar
%  disp(sprintf('addcb(''hor'',%s,%s,%s)',...
%	       inputname(5),inputname(5),inputname(3))) 
end

% Prepare output
varns={h,cax,X,c11,cmn};
varargout=varns(1:nargout);
