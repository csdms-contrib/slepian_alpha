function varargout=fillbox(lrtbmat,colvec,varargin)
% plothandles=FILLBOX(lrtbmat,colvec,opt)
%
% lrtbmat   mX4 for m boxes whose Left Right Tob Bottom you give
% colvec    an mX3 RGB matrix
%           an mX1 or 1Xm vector of strings
%           1 scalar (automatic indexing) or 1 letter string
%           an mX1 or 1Xm vector of scalars
% opt       0 planar plot [default]
%           1 makes a three-dimensional plot
%
% EXAMPLE I:
%
% lrtb=[50 100 80 100];
% load clown; image(X); hold on; a=fillbox(lrtb,'w')
% b=boxmid(lrtb); c=text(b(1),b(2),'Clown');
% set(c,'HorizontalA','center','FontS',15)
% set(a,'FaceColor','y','EdgeColor','b')
%
% EXAMPLE II:
%
% lrtbmat=[10 20 20 10 ; 20 30 21 11 ; 31 41 29 22 ; 12 14 3 0];
% fillbox(lrtbmat,2)
% fillbox(lrtbmat,[1 2 3 4])
% fillbox(lrtbmat,'y')
% fillbox(lrtbmat,'ykbr')
% fillbox(lrtbmat,['y' 'k' 'b' 'r'])
% fillbox(lrtbmat,[{'y'},{'k'},{'b'},{'r'}])
% fillbox(lrtbmat,rand(size(lrtbmat,1),3))
%
% FILL3 gets fed MXN matrix for N polygons of order M
%
% See also: EXT2LRTB
%
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2012

defval('colvec',grey(8))

if size(lrtbmat,2)~=4; error('LRTB not mX4'); end

if ischar(colvec) | iscell(colvec)
  colvec=getcol(colvec);
  if size(colvec,1)==1; colvec=repmat(colvec,size(lrtbmat,1),1);end 
else
  % If it's a matrix
  if prod(size(colvec))~=1 % Not just a number
    if (size(colvec,2)~=3 | size(colvec,1)~=size(lrtbmat,1)) ...
	  & prod(size(colvec))~=size(lrtbmat,1)
      error('Input color matrix is not in right format')
    end    
  end
end

% Need to get it to give you a color specification.
if prod(size(colvec))~=size(lrtbmat,1)
  colvec=colvec';
else
  colvec=colvec(:)';
end

[m,n]=size(lrtbmat);
[left,right,top,bottom]=...
    deal(lrtbmat(:,1),lrtbmat(:,2),lrtbmat(:,3),lrtbmat(:,4));

EKS=[left left right right]';
WAI=[top bottom bottom top]';

indo=repmat(1:size(EKS,2)',3,1);
indo=indo(:)';

% Plotting
if nargin>2 
  if prod(size(colvec))==1 | size(colvec,1)==1
    plothandles=fill3(EKS,WAI,repmat(1,4,size(EKS,2)),colvec);
  else
    stev=sprintf('EKS(:,%i),WAI(:,%i),[1 ; 1; 1; 1],colvec(:,%i)'',',indo);
    eval(['plothandles=fill3(' stev(1:end-1) ');'])
  end
else   
  if prod(size(colvec))==1| size(colvec,1)==1
    plothandles=fill(EKS,WAI,colvec);
  else
    stev=sprintf('EKS(:,%i),WAI(:,%i),colvec(:,%i)'',',indo);
    eval(['plothandles=fill(' stev(1:end-1) ');'])
  end
end

% Produce output
varns={plothandles};
varargout=varns(1:nargout);
