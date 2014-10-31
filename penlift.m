function varargout=penlift(varargin)
% [X,Y,Z,p]=penlift(X,Y,Z,dlev)
% [X,Y]=penlift(X,Y,dlev)
% XYZ=penlift(XYZ,dlev)
%
% Lifts the pen by inserting NaNs where the jumps are deemed to big
%
% INPUT:
%
% X,Y,Z      Coordinates (one Mx3 matrix may replace three Mx1 matrices)
% dlev       Whatever exceeds the metric du jour [default: 3]
%
% OUTPUT:
%
% X,Y,Z      New coordinate matrices, format of input matched
% p          The special points at which the pen was being lifted
%
% Last modified by fjsimons-at-alum.mit.edu, 01/26/2012

if nargin==1 || [nargin==2 && prod(size(varargin{2}))==1]
  % Input is (XYZ) or (XYZ,dlev)
  XYZ=varargin{1};
  X=XYZ(:,1);
  Y=XYZ(:,2);
  flag=0;
  if size(XYZ,2)==3
    Z=XYZ(:,3);
    clear XYZ
    flag=1;
  end
  if nargin==2
    % Input is indeed (XYZ,dlev)
    dlev=varargin{2};
  end
elseif nargin>=2
  % Input is (X,Y,...) or (X,Y,Z,...)
  X=varargin{1};  
  Y=varargin{2};
  difer(size(X)-size(Y),[],[],NaN)
  flag=2;
  % Input is (X,Y,Z) or (X,Y,dlev)
  if nargin>=3 && prod(size(varargin{3}))~=1
    Z=varargin{3};
    difer(size(X)-size(Z),[],[],NaN)
    flag=3;
    % input is (X,Y,Z,dlev)
    if nargin==4
      dlev=varargin{4};
    end
  elseif nargin==3
    % input is indeed (X,Y,dlev)
    dlev=varargin{3};
  end
end

% Still may have had an empty level, and there may be no Z
defval('dlev',3)
defval('Z',[])

% Distance between two consecutive points in the plane
d=sqrt([diff(X)].^2+[diff(Y)].^2);

% Now right after each of these positions need to insert a NaN
if prod(size(X))==length(X)
  % They're all simply vectors

  % Where do the NaN's go (linear index is enough here)
  p=find(d>dlev*nanmean(d));

  X=insert(X,NaN,p+1);
  Y=insert(Y,NaN,p+1);
  if ~isempty(Z)
    Z=insert(Z,NaN,p+1);
  end
else
  % They're all same-sized matrices - but they may not require NaN's in
  % all the same places. If they don't, must return them as cells, which
  % we are going to assume to be sure; also, we make them bigger
  
  % Where do the NaN's go (row and column index needed)
  [p,j]=find(d>dlev*repmat(nanmean(d),size(d,1),1));

  % Every column is a different curve so to speak
  for in=1:size(X,2)
    XX{in}=insert(X(:,in),NaN,p(j==in)+1);
    YY{in}=insert(Y(:,in),NaN,p(j==in)+1);
    if ~isempty(Z)
      ZZ{in}=insert(Z(:,j(in)),NaN,p(j==1)+1);
    end
  end
  % Now reassign the right variable name
  X=XX; Y=YY; if ~isempty(Z); Z=ZZ; end
end

% Output preparation
if flag==0 || flag==1
  XYZ(:,1)=X;
  XYZ(:,2)=Y;
  if flag==1
    XYZ(:,3)=Z;
  end
  vars={XYZ};
else
  vars={X(:),Y(:),Z(:),p};
end

% Actual output generation
varargout=vars(1:nargout);
