function RGB=joincolmap(A,RGB,I,colmap,cax,nancol,ncol)
% RGB=joincolmap(A,RGB,I,colmap,cax,nancol,ncol)
%
% Data matrix 'A' of size MxN
% Matrix 'RGB' of size MxNx3 colors
% Index matrix 'I' or matrix MxN of logicals
% Colormap 'colmap' with limits 'cax' and a 3x1 matrix
% with the color NaN data values will get
% Scalar 'ncol' is the resolution of the colormap
%
% From karason-at-alum.mit.edu

if [isempty(I) & isempty(RGB)] | ...
      [~isempty(I) & isempty(RGB)] | ...
      [~isempty(I) & ~isempty(RGB)] 
  % Check if 'ncol' exists as a variable
  if exist('ncol','var'),
    cont=1;
  else
    cont=0;
  end
  [n,m]=size(A);
  
  if islogical(I)
    I=find(I);
  end
  
  defval('I',(1:n*m)');
  defval('RGB',repmat(nan,[size(A) 3]));
  defval('ncol',linspace(cax(1),cax(2),size(colmap,1))');
  defval('nancol',[1 1 1]);
  
  A=A(I(:));
  
  thenans=find(isnan(A));
  A(thenans)=0;
  A(A<cax(1))=cax(1);
  A(A>cax(2))=cax(2);
  
  if cont
    ncol(ncol<cax(1))=cax(1);
    ncol(ncol>cax(2))=cax(2);
    %cols=interp1(linspace(cax(1),cax(2),size(colmap,1)),colmap,ncol(:),'nearest');
    cols=interp1(ncol,colmap,ncol(:),'nearest');
    koler=repmat(nan,[length(A) 3]);
    for j=1:size(cols,1),
      Il=find(A>=ncol(j));
      if ~isempty(Il),
	koler(Il,:)=repmat(cols(j,:),size(Il));
      end
    end  
  else
    koler=interp1(ncol,colmap,A,'nearest');
  end
    
  if ~isempty(thenans),
    koler(thenans(:),1)=nancol(1);
    koler(thenans(:),2)=nancol(2);
    koler(thenans(:),3)=nancol(3);
  end
   
  RGB(I)=koler(:,1);
  RGB(I+n*m)=koler(:,2);
  RGB(I+2*n*m)=koler(:,3);
else  
  if isempty(I)
    RGB=RGB;
  end  
end
