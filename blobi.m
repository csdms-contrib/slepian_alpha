function varargout=blobi(Nj,plotit)
% varargout=blobi(Nj,plotit)
%
% Makes an embedded indicator function for a random blob.
%
% INPUT:
%
% Nj       Smoothness measure [default: 3]
% plotit   1 Will make a plot
%          0 Will not make a plot
%
% OUTPUT:
%
% I           The indicator function 
%
% Last modified by fjsimons-at-alum.mit.edu, 01/26/2023


% Defaults
defval('Nj',3)
defval('plotit',0)

% Random smoothish blob, only make one, really could use RANDCIRC, but hey
[xb,yb]=blob(1,Nj);
% Embedding and grid coarseness
emb=4; grd=1;
rx=range(xb); ry=range(yb);
% Get some idea of the blob-dependent thus random spacing
%x=linspace(min(xb)-rx/4,max(xb)+rx/emb,length(xb)/2);
%y=linspace(min(yb)-ry/4,max(yb)+ry/emb,length(yb)/2);
% The above two lines always yield the same dimension
x=min(xb)-rx/4:grd*median(abs(diff(xb))):max(xb)+rx/emb;
y=min(yb)-ry/4:grd*median(abs(diff(yb))):max(yb)+ry/emb;
if isempty(x) || isempty(y); kb; error('refine grid'); end
% Make grid
[X,Y]=meshgrid(x,y);
% Find the inside of the blob
I=inpolygon(X,Y,xb,yb);

if plotit==1
    imagesc(x,y,I)
    hold on
    plot(xb,yb,'k','LineWidth',2)
    hold off
    axis image xy
    longticks(gca,2)
    t=title(sprintf('%i x %i / %i%% coverage',...
                    size(I,1),size(I,2),round(100*sum(I(:))/prod(size(I)))));
    movev(t,ry/emb/10)
end

% Optional output
varns={I};
varargout=varns(1:nargout);
