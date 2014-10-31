function varargout=eqband(res,inclang)
% XY=eqband(res,inclang)
% eqband(...) % Only makes a plot
%
% Finds the coordinates of a band on the sphere
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% inclang  Inclination angle, determines size polar caps
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
%
% SEE ALSO: ALLOCEANS
%
% Written by D.C. Slobbe, 21/10/2010
% Last modified by fjsimons-at-alum.mit.edu, 11/1/2010

defval('res',0)
defval('inclang',66.15);

% The directory where you keep the coordinates
whereitsat=fullfile(getenv('IFILES'),'COASTS');

if res==0
  fnpl=fullfile(whereitsat,sprintf('eqband-%5.2f.mat',inclang));
else
  fnpl=fullfile(whereitsat,sprintf('eqband-%5.2f-%i.mat',inclang,res));
end

if exist(fnpl,'file')==2
  load(fnpl)
  if nargout==0
    plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
  else
    varns={XY};
    varargout=varns(1:nargout);
  end
else
  if res==0
    % Resolution seems overkill to me, but let's leave it
    NX=3599;
    X=linspace(0.1,359.9,NX)';
    NY=round(inclang/90*NX/2)+1;
    Y=linspace(inclang,-inclang,NY)';
    XY=[[X ; ones(NY,1)*360 ; flipud(X) ; zeros(NY,1)            ],...
	[ones(NX,1)*-inclang ; flipud(Y) ; ones(NX,1)*inclang ; Y]];
    XY(:,1)=XY(:,1)+360;
    
    % Eyeball
    plot(XY(:,1),XY(:,2),'LineW',2,'Color','k');
    axis equal
    
    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curvecheck(XY(:,1),XY(:,2)) 
    
    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(fnpl,'XY')
  else
    XY =eqband(0);
    XYb=bezier(XY,res);
    XY=XYb;
    save(fnpl,'XY')
  end
  varns={XY};
  varargout=varns(1:nargout);
end
