function varargout=plotonchunk(xi,eta)
% [p,pg]=PLOTONCHUNK(xi,eta)
% 
% Plots scattered data on the faces of the cubed sphere in the plane
%
% INPUT:
%
% xi,eta   6-dimensional cell arrays with points to be plotted
%          In the two-dimensional rendition, xi goes up and corresponds
%          to the original longitude of the x+ chunk, while eta goes
%          across to the right and corresponds to the colatitude of the
%          x+ chunk 
%
% OUTPUT:
%
% p,pg     Axis handles to the plot lines and the grid
%
% SEE ALSO: PLOTONCUBE, PLM2CUBE (notably, for demos using this)
%
% Last modified by fjsimons-at-alum.mit.edu, 08/19/2009

for in=1:6
  xup=floor(in/2)*pi/2;
  yup=(ceil(in/2)-1)*pi/2;
  if ~isempty(eta{in})
    p(in)=plot(eta{in}+xup,xi{in}+yup,'.');
  else
    p(in)=NaN;
  end
  hold on
  pg{in}=fridplot([xup-pi/4 pi/4+xup],[yup-pi/4 pi/4+yup]);
end

hold off
set(p(~isnan(p)),'MarkerS',2)
axis off
%set(gcf,'color','w')

% Prepare output
varns={p,pg};
varargout=varns(1:nargout);

