function varargout=blob(N,Nj)
% [x,y]=BLOB(N,nj)
%
% Makes (moving) picture of smoothly deforming random blobs. 
%
% INPUT:
%
% N     Number of loops for, and if movie [default: 100]
% Nj    Smoothness, roughly [default: 10]
%
% OUTPUT:
%
% x,y    Horizontal and vertical coordinates. Figure is being created
%        only no output requested.
%
% SEE ALSO:
%
% RANDCIRC
%
% Tested on 8.3.0.532 (R2014a), 9.0.0.341360 (R2016a), 9.7.0.1190202 (R2019b)
%
% Last modified by fjsimons-at-alum.mit.edu, 01/26/2023

% Defaulted
defval('N',100)
defval('Nj',10);
% Makes the ur-random circle
[xold,yold]=randcirc(0,0,1,1,10);
% Discretize the radius perturbation
r=linspace(0,1,Nj+1);
rm=1-r;
if nargout==0
  figure(gcf)
end

% Loop and plot
for index = 1:N
  [x,y]=randcirc(0,0,1,0.2,10);                                             
  col=rand(1,3);
  for j=1:Nj,
    % First is increasing radius, second decreasing
    xx=x*r(j)+xold*rm(j);
    yy=y*r(j)+yold*rm(j);      
    if nargout==0
%      plot(xx,yy,'-','linewidth',2,'color',[1,1,1]*0);
%      plot(xx,yy,'-','linewidth',2,'color',col);
      axis([-2 2 -2 2])
      axis off
      drawnow
%      pause(0.5)
      fill(xx',yy',col)
      axis equal ; axis off   
    end
  end
  xold=x;
  yold=y;
end

% Optional output
varns={xx(:),yy(:)};
varargout=varns(1:nargout);
  
