function revolved=revolve(radfun,meth1,meth2,init)
% revolved=REVOLVE(radfun,meth1,meth2,init)
%
% Makes a matrix with the function values in 'radfun' as radial profile 
% in any direction trough the center. Interpolation methods used are 'meth1' for
% initial interpolation, 'meth2' for later interpolations.
% Defaults are 'spline' and 'linear'. Function 'radfun' must represent
% equally spaced data. 
% The initial resolution is given in 'init', for which the default is 5.
% The algorithm first interpolates 'radfun' 'init' times with 'meth1'
% and then uses a simples interpolation for every radial profile through the
% mxm matrix that will be created, if m=length(radfun).
%
% EXAMPLE:
%
%     m=50; nlev=0.4;
%     AA=(1-nlev*rand(1,m)).*sin(linspace(pi/2,pi,m));
%     subplot(121) ; plot(AA)
%     revolved=revolve(AA);
%     subplot(122) ; surf(revolved) ; shading flat
%
% Last modified by fjsimons-at-alum.mit.edu, May 16th 2001

if nargin==1
  disp('Using default interpolation schemes')
  meth1= 'spline';
  meth2= 'linear';
  init=5;
end

radfun=radfun(:);
m=length(radfun);
lindom=linspace(1,m,m*init);
radint=interp1(1:m,radfun,lindom,meth1);

% Create 1/4 of the matrix domain
[X,Y]=meshgrid(1:m,1:m);

% Distance from origin (0,0)
dist=sqrt((X-1).^2+(Y-1).^2)+1;

funint=reshape(interp1(lindom,radint,dist(:),meth2),m,m);

revolved=flipflop(fliplr(flipud(funint)));
