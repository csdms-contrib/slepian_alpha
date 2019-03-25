function varargout=alloceans3(res,inclang)
% XY=ALLOCEANS3(res)
% ALLOCEANS(...) % Only makes a plot
%
% Finds the coordinates of some of the worlds' oceans so we can combine
% its localization kernel with those for the missing continents to turn
% into Slepian eigenfunctions for all of the world's oceans
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
% Last modified by fjsimons-at-alum.mit.edu, 11/1/2010
% Last modified by charig-at-princeton.edu, 08/09/2015

res=10;
defval('inclang',50);

if inclang==90
  fnpl=fullfile(getenv('IFILES'),'COASTS',sprintf('alloceans3-%i.mat',res));
else
  fnpl=fullfile(getenv('IFILES'),'COASTS',sprintf('alloceans3-%i-%.2f.mat',res,inclang));
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
  
    ShftAngl    = 180;
    %XYbos       = eqband(0,inclang);
    %XYbos(:,1)  = XYbos(:,1)-360-ShftAngl;

    XY_NA       = namerica(res);
    XY_NA(:,1)  = XY_NA(:,1);

    XY_GR       = greenland(res,[],0);
    XY_GR(:,1)  = XY_GR(:,1);

    XY_EA       = eurasia(res);
    XY_EA(:,1)  = XY_EA(:,1);

    XY_AF       = africa(res);
    XY_AF(:,1)  = XY_AF(:,1);
    
    XY_SA       = samerica(res);
    XY_SA(:,1)  = XY_SA(:,1);
    
    XY_AU       = australia(res);
    XY_AU(:,1)  = XY_AU(:,1);
    
%     % Antarctica receives special treatment
%     clf
%     [~,handl,XYant]=plotcont([0 -62],[360 -83],[],0);
%     delete(handl)
% 
%     % Get rid of common NaNs
%     XYant      = XYant(~isnan(XYant(:,1)) & ~isnan(XYant(:,2)),:);
%     XYant=wrapTo180(XYant);
%     XYant(1,1)=-180;
%     XYant=XYant([1:227 230:end],:);
    XYant=antarctica(res,[],1);

    % Form closed contour
    %XYant      = XYant(~isnan(XYant(:,1)) & ~isnan(XYant(:,2)),:);
    %XYant      = XYant([231:end,1:227],:);
    
    [XY_NA(:,1),XY_NA(:,2)] = poly2cw(XY_NA(:,1), XY_NA(:,2));
    [XY_GR(:,1),XY_GR(:,2)] = poly2cw(XY_GR(:,1), XY_GR(:,2));
    [XY_EA(:,1),XY_EA(:,2)] = poly2cw(XY_EA(:,1), XY_EA(:,2));
    [XY_AF(:,1),XY_AF(:,2)] = poly2cw(XY_AF(:,1), XY_AF(:,2));
    [XYant(:,1),XYant(:,2)] = poly2cw(XYant(:,1), XYant(:,2));
    [XY_SA(:,1),XY_SA(:,2)] = poly2cw(XY_SA(:,1), XY_SA(:,2));
    [XY_AU(:,1),XY_AU(:,2)] = poly2cw(XY_AU(:,1), XY_AU(:,2));

    % Now combine these into an overall interior coastline
    %warning off map:vectorsToGPC:noExternalContours
    XY=[XY_NA; NaN NaN; XY_GR; NaN NaN; XY_EA; NaN NaN; XY_AF; NaN NaN; XY_SA; NaN NaN; XY_AU; NaN NaN; flipud(XYant)];
    %XY=wrapTo180(XY);
    [tempy,tempx]=flatearthpoly(flipud(XY(:,2)),flipud(XY(:,1)),180);
    XY=[tempx tempy];
    %warning on map:vectorsToGPC:noExternalContours

    % This should be done manually
    %StrtIndpol  = [0;find(isnan(XY(:,1)));length(XY(:,1))+1];
    %[~,maxdi]   = max(diff(StrtIndpol));
    %XY          = [XY(StrtIndpol(maxdi)+1:StrtIndpol(maxdi+1)-1,:)];
    
    
    %XY(:,1)     = XY(:,1)+360;

    % Eyeball
    plot(XY(:,1),XY(:,2),'LineW',2,'Color','k');
    axis equal
    
    
    if inclang ~= 90
        test=polarcaps(0,inclang);
        [x1,y1] = polybool('subtraction',XY(:,1),XY(:,2),test(:,1),test(:,2));
        XY = [x1 y1];
    end
    
    [LATo,LONo] = interpm(XY(:,2),XY(:,1),0.25);
    XY          = [LONo,LATo];
   
    
    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %curvecheck(XY(:,1),XY(:,2)) 

    % Check this out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(fnpl,'XY')

  varns={XY};
  varargout=varns(1:nargout);
end
