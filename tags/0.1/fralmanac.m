function val=fralmanac(neem,plenet)
% val=FRALMANAC('name','plenet')
%
% Accesses a database with planetary constants given in SI units.
% If 'name' is empty returns a list of possibilities.
%
% INPUT: 
%
% 'name'    -> for plenet 'Earth' [default]
%
%           'CMB'           Radius of the core-mantle-boundary  [m]
%           'DegDis'        Length of equatorial longitude      [m]
%           'GM_EGM2008'    EGM2008 reference mass constant     [m^3s^{-2}]
%           'GM_EGM96'      EGM96 reference mass constant       [m^3s^{-2}]
%           'GM_EIGENCG03C' EIGENCG03C reference mass constant  [m^3s^{-2}]
%           'GravAcc'       Surface gravity                     [ms^{-2}]
%           'GravCst'       Gravitational constant              [m^3kg^{-1}s^{-2}]
%           'IMR2'          Reduced moment of inertia (I/MR2)   [dimensionless]
%           'Mass'          Reference mass                      [kg]
%           'Radius'        Volumetric mean radius              [m]
%           'a_EGM2008'     EGM2008 reference radius            [m]
%           'a_EGM96'       EGM96 reference radius              [m]
%           'a_EIGENCG03C'  EIGENCG03C reference radius         [m]
%           'omega_wgs84'   Rotational velocity                 [s^-1]
%           'rf_wgs84'      Inverse flattening                  [s^-1]
%           'a_wgs84'       Semimajor axis                      [m]
%           'GM_wgs84'      WGS84 reference mass constant       [m^3s^{-2}]
%
%           -> for plenet 'Mars'
% 
%           'DegDis'        Length of equatorial longitude      [m]
%           'GM_GMM2B'      GMM2B reference mass constant       [m^3s^{-2}]
%           'GM_JGM85H02'   JGM85H02 reference mass constant    [m^3s^{-2}]
%           'GravAcc'       Surface gravity                     [ms^{-2}]
%           'IMR2'          Reduced moment of inertia (I/MR2)   [dimensionless]
%           'Mass'          Reference mass                      [kg]
%           'Radius'        Volumetric mean radius              [m]
%           'a_GMM2B'       GMM2B reference radius              [m]
%           'a_JGM85H02'    JGM85H02 reference radius           [m]
%
%           -> for plenet 'Moon'
%
%           'DegDis'        Length of equatorial longitude      [m]
%           'GM_GLGM2'      GLGM2 reference mass constant       [m^3s^{-2}]
%           'GravAcc'       Surface gravity                     [ms^{-2}]
%           'IMR2'          Reduced moment of inertia (I/MR2)   [dimensionless]
%           'Mass'          Reference mass                      [kg]
%           'Radius'        Volumetric mean radius              [m]
%           'a_GLGM2'       GLGM2 reference radius              [m]
%
%           -> for plenet 'Venus'
%
%           'a_SHTJV360'    SHTJV360 v A02 reference radius     [m]
%           'a_SHGJ180U'    SHGJ180U v A01 reference radius     [m]
%
%           -> for plenet 'SHM', i.e. spherical harmonics coefficients
% 
%           'EGM2008_Topography' Earth topography (to degree 2190)
%           'EGM2008_ZeroTide'   Earth geopotential (to degree 2190)
%           'EGM96'              Earth geopotential (to degree 360)
%           'EIGEN_CG03C'        Earth geopotential (to degree 360)
%           'GLGM2'              Moon geopotential (from Clementine)
%           'GLTM2B'             Moon topography (from Clementine)
%           'GMM2B'              Mars geopotential (from MGS / Frank Lemoine)
%           'GTM090'             Mars shape (from MOLA / Greg Neumann)
%           'GTM3AR'             Earth topography (from Georg Wenzel)
%           'JGM85H02'           Mars geopotential (from MGS)
%           'Mars2000'           Mars topography (from MOLA / Greg Neumann)
%           'MarsTopo719'        Mars shape (from Mark Wieczorek)
%           'ULCN359_lpo'        Lunar shape (from Mark Wieczorek)
%           'VenusTopo719'       Venus shape (from Mark Wieczorek)
%           'SHTJV360'           Venus topo/radius, version A02 (via Kevin Lewis)
%           'SHGJ180U'           Venus gravitational potential (via Kevin Lewis)
% 
%           -> for plenet 'XYZ', i.e. globally expanded planetary models
%
%           'EGM96'              Earth free-air gravity (to degree 360)
%           'GLGM2not02'         Moon free-air gravity, no degree 2 (from Clementine)
%           'GLTM2B'             Moon topography (from Clementine)
%           'GTM090'             Mars shape (from MOLA / Greg Neumann)
%           'GTM3AR'             Earth topography (from Georg Wenzel)
%           'Mars2000'           Mars topography (from MOLA / Greg Neumann)
%           'MarsTopo719'        Mars shape (from Mark Wieczorek)
%           'MarsTopo719not02'   Mars shape, no degree 2
%           'ULCN359_lpo'        Lunar shape (from Mark Wieczorek)
%           'ULCN359_lponot02'   Lunar shape, no degree 2
%           'VenusTopo719'       Venus shape (from Mark Wieczorek)
%           'VenusTopo719not02'  Venus shape, no degree 2
%
% 'plenet'  'Earth' [default]
%           'Mars'
%           'Moon'
%            ...the above contain the basic constants in the first list
%           'SHM' spherical harmonics models
%           'XYZ' spatially expanded models
%            ...the above contain the fields in the second list
%
% EXAMPLE:
%
% fralmanac([],'Moon') % Tells you what you can get for this planet
% fralmanac([],'SHM') % Tells you what you can get for this class of models
% v=fralmanac('Radius','Mars'); % Returns a radius for Mars
% v=fralmanac('EGM96','XYZ'); % Returns free-air gravity for Earth
% v=fralmanac('EGM96','SHM'); % Returns geopotential coefficients for Earth
% v=fralmanac('SHTJV360','SHM'); % Returns Venus topography
%
% SEE ALSO:
%
% ORDERFIELDS
%
% TO DO:
% 
% Should rename the expanded models more specifically than is the case now.
%
% Among others, from http://nssdc.gsfc.nasa.gov/planetary/factsheet/
%
% Last modified by fjsimons-at-alum.mit.edu, 02/22/2012

defval('neem',[])
defval('plenet','Earth')

% Change this to load the variable instead of the whole file
load(fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS',plenet))

if nargin>0 & ~isempty(neem) 
    val=eval([plenet '.' neem]);
    if iscell(val) & prod(size(val))==1
      val=val{1};
    end
else
  number(str2mat(fieldnames(eval(plenet))))
end
