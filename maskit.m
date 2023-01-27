function varargout=maskit(v,p,scl)
% [v,cr,I]=MASKIT(v,p,scl)
%
%
% INPUT:
%
% v       A vector that is the unwrapping of a matrix
% p       A parameter structure with AT LEAST this constant
%           NyNx  number of samples in the y and x directions
%           mask  an index matrix with a mask, OR:
%                 a region name that will be scaled to fit
% scl     A scale between 0 and 1 with the occupying center fraction;
%         due to strict inequality at maximum always leaves one pixel rim
%
% OUTPUT:
%
% v       The masked output matrix, unrolled into a vector
% cr      The colum,row index coordinates of the masking curve
% I       The mask, unrolled into a vector
%
% EXAMPLE:
%
% p.mask='france';
% p.quart=0; p.blurs=Inf; p.kiso=NaN; clc; 
% [Hx,th,p]=simulosl([],p,1); 
% [Hm,cr,I]=maskit(Hx,p);
% subplot(121); imagefnan([1 1],p.NyNx([2 1]),v2s(Hx,p),[],halverange(Hx,80)); axis ij image; 
% hold on; twoplot(cr,'Color','k'); hold off; longticks; grid on
% xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)]) 
% subplot(122); imagefnan([1 1],p.NyNx([2 1]),v2s(Hm,p),[],halverange(Hx,80)); axis ij image
% hold on; twoplot(cr,'Color','k'); hold off; longticks; grid on
% xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)])
% figdisp(sprintf('maskit_%s',p.mask),[],[],2)
%
% Last modified by fjsimons-at-alum.mit.edu, 01/26/2023

% Get the mask
if isstr(p.mask)
    XY=eval(p.mask);
end

% Default scale is 80% of the area
defval('scl',0.8)

% Scale the mask, remember for the curve the Y goes up, for the image, down
cr(:,1)=          scale(XY(:,1),[1 p.NyNx(2)]+[1 -1]*(1-scl)/2*(p.NyNx(2)-1));
cr(:,2)=p.NyNx(1)-scale(XY(:,2),[1 p.NyNx(1)]+[1 -1]*(1-scl)/2*(p.NyNx(1)-1))+1;

% The underlying image grid
[X,Y]=meshgrid(1:p.NyNx(2),1:p.NyNx(1));

% Determine the mask
I=inpolygon(X,Y,cr(:,1),cr(:,2));

% Apply the mask
v(~I(:))=NaN;

% Variable output
varns={v,cr,I};
varargout=varns(1:nargout);
