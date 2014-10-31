function varargout=alloceans(res)
% XY=ALLOCEANS(res)
% ALLOCEANS(...) % Only makes a plot
%
% Finds the coordinates of some of the worlds' oceans so we can combine
% its localization kernel with those for the missing continents to turn
% into Slepian eigenfunctions for all of the world's oceans. Also, there
% is a polar gap taken out. 
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% inclang  Inclination angle, determines size of polar caps
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the world's ocean
%
% SEE ALSO: EQBAND
%
% Written by D.C. Slobbe, 21/10/2010
% Last modified by fjsimons-at-alum.mit.edu, 10/08/2011

defval('res',0)
defval('inclang',66.15);

if res==0
  fnpl=fullfile(getenv('IFILES'),'COASTS','alloceans.mat');
else
  fnpl=fullfile(getenv('IFILES'),'COASTS',sprintf('alloceans-%i.mat',res));
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
    ShftAngl    = 180;
    XYbos       = eqband(0,inclang);
    XYbos(:,1)  = XYbos(:,1)-360-ShftAngl;

    XY_NA       = namerica(0);
    XY_NA(:,1)  = XY_NA(:,1)-360;

    XY_GR       = greenland(0);
    XY_GR(:,1)  = XY_GR(:,1)-360;

    XY_EA       = eurasia(0);
    XY_EA(:,1)  = XY_EA(:,1)-360;

    % Antarctica receives special treatment
    clf
    [~,handl,XYant]=plotcont([0 -62],[360 -83],[],0);
    delete(handl)

    % Get rid of common NaNs
    XYant      = XYant(~isnan(XYant(:,1)) & ~isnan(XYant(:,2)),:);

    % Form closed contour
    XYant      = XYant(~isnan(XYant(:,1)) & ~isnan(XYant(:,2)),:);
    XYant      = [[0,-90];[0,XYant(227,2)];...
		  XYant([231:end,1:227],:);[360,-90];[0,-90]];
    [XY_NA(:,1),XY_NA(:,2)] = poly2cw(XY_NA(:,1), XY_NA(:,2));
    [XY_GR(:,1),XY_GR(:,2)] = poly2cw(XY_GR(:,1), XY_GR(:,2));
    [XY_EA(:,1),XY_EA(:,2)] = poly2cw(XY_EA(:,1), XY_EA(:,2));
    [XYant(:,1),XYant(:,2)] = poly2cw(XYant(:,1), XYant(:,2));

    % Now combine these into an overall interior coastline
    warning off map:vectorsToGPC:noExternalContours
    [Xu,Yu]=polybool('subtraction',XYbos(:,1),XYbos(:,2),...
		     XY_NA(:,1),XY_NA(:,2)); ...
    [Xu,Yu]=polybool('subtraction',Xu,Yu,XY_GR(:,1),XY_GR(:,2));
    [Xu,Yu]=polybool('subtraction',Xu,Yu,XYant(:,1),XYant(:,2));
    [Xu,Yu]=polybool('subtraction',Xu,Yu,XY_EA(:,1),XY_EA(:,2));
    % Make sure you get both pieces
    XY_EA(:,1)=XY_EA(:,1)-360;
    [Xu,Yu]= polybool('subtraction',Xu,Yu,XY_EA(:,1),XY_EA(:,2));
    XYant(:,1)=XYant(:,1)-360;
    [Xu,Yu]=polybool('subtraction',Xu,Yu,XYant(:,1),XYant(:,2));
    XY=[Xu,Yu];
    warning on map:vectorsToGPC:noExternalContours

    % This should be done manually
    StrtIndpol  = [0;find(isnan(XY(:,1)));length(XY(:,1))+1];
    [~,maxdi]   = max(diff(StrtIndpol));
    XY          = [XY(StrtIndpol(maxdi)+1:StrtIndpol(maxdi+1)-1,:)];
    [LATo,LONo] = interpm(XY(:,2),XY(:,1),0.25);
    XY          = [LONo,LATo];
    
    XY(:,1)     = XY(:,1)+360;

    % Eyeball
    plot(XY(:,1),XY(:,2),'LineW',2,'Color','k');
    axis equal
    
    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curvecheck(XY(:,1),XY(:,2)) 

    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(fnpl,'XY')
  else
    XY=alloceans(0);
    XYb=bezier(XY,res);
    XY=XYb;
    save(fnpl,'XY')
  end
  varns={XY};
  varargout=varns(1:nargout);
end
