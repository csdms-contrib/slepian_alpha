function varargout = potupxyz(coefin,rold,lon,lat,rnew,Lmax,inv,Lrange)
% dataout=potupxyz(coefin,rold,lon,lat,rnew,Lmax,inv,Lrange)
%
% Upward-continues potential-field scalar spherical harmonic coefficients
% to a series of xyz (lon,lat,radius) locations. This can simulate in situ
% satellite observations.
% 
%
% INPUT:
%
% coefin    The field (at rold) you wish to evaluate:
%           [1] Spherical harmonic coefficient vector, OR 
%           [2] Matrix whose columns are coefficients for different fields 
%           you wish to evaluate, OR
%           [3] A single field in lmcosi format.
%            -In the case of 1 and 2, the lm cycling must be in the columns
%             and the vectors must have orders -l to l as in addmout
%             because we use ylm. Lmcosi for option 3 must be 
%             complete 0 to Lmax.
% rold      The current radius for your coefficients. [default: Earth 
%             radius for coeficients of potential field]
% lon       "lon": a column vector with longitudes [degrees]
% lat       "lat": a column vector with latitudes [degrees]
% rnew      A column vector of new radii that correspond to lon and 
%             lat vectors. [meters]
% Lmax      maximum degree [no default]
% inv       up (0) or down (1)? (Handbook: A->0 , A^{-1}->1 )
%           up also means derivative, down means integral [default: 0]
% Lrange    all the L that are used if not all l,m combinations between 0 
%           and Lmax are included in the matrix (for example in sdwcapup)
%           WARNING: Use correct Ls. Otherwise you seriously mess things
%           up.
%
% OUTPUT:
%
% dataout   Values of the initial field evaluated at the coordinates. One
%           column for each field you gave.
%
% Last modified by charig-at-arizona.edu, 06/27/2022
% Last modified by plattner-at-alumni.ethz.ch, 05/09/2018

%defval('coefin',fralmanac('EGM2008_Topography','SHM'))
defval('rold',fralmanac('a_EGM96'))
defval('rnew',(rold+1))
defval('Lrange',[])
defval('inv',0)
% There should be no defval for Lmax


if ~ischar(coefin) && (addmoff(Lmax) == size(coefin,1) || ~isempty(Lrange))
    % We have column vectors of spherical harmonic coefficients that start
    % at zero, OR we have a Lrange (which may or may not be continuous)

    % How many points are we evaluating?
    n = length(lon);

    if ~isempty(Lrange)
        %disp('POTUPXYZ: Using provided L range')
        bigl=Lrange(:);
    else
        [~,bigl] = addmout(Lmax);
    end

    % Get the spherical harmonic degrees and orders
    %[EM,EL]=addmout(L);

    % Define an upward continuation operator for the 
    % spherical harmonics
    bigl_a=repmat(bigl,1,n);
    A = (repmat((rnew/rold)',length(bigl),1)).^(-bigl_a - 1);
    
    % Get The spherical-harmonics sensitivity matricies for our orbital points
    % And change them from unit normalized to 4pi normalized harmonics
    [YA,~,~,ems,ells]=ylm([0 Lmax],[],pi/2-lat*pi/180,lon*pi/180,[],[],[],1);
    YA=YA*sqrt(4*pi);
    % Don't forget the ems to correct the phase
    YA = YA.*repmat((-1).^ems,1,n);

    % Since YLM can only do all of the degrees not subsets, we now remove
    % the extra ones if we had a Lrange
    indks = any(ells==unique(bigl'),2);
    YA = YA(indks,:);

    % Now multiply the SH sensitivities by the continuation operator and by
    % the field coefficients, to get the new fields evluated at our points
    % Don't forget the ems
    if (size(coefin,1)==length(bigl))
        for i=1:size(coefin,2) 
            dataout(:,i) = sum(A.*YA.*coefin(:,i),1)';
        end
    else
        error('Number of coefin entries does not match Lrange')
    end

    varns={dataout};
    varargout=varns(1:nargout);


elseif ~ischar(coefin) && (addmup(Lmax) == size(coefin,1))
    % We have a single lmcosi matrix, which is not missing any degrees.
    % Reorder it into a column and use the code above.
    [~,~,~,~,~,~,~,~,~,ronm]=addmon(Lmax);
    temp = coefin(:,3:4);
    coefin = temp(ronm);

    dataout = potupxyz(coefin,rold,lon,lat,rnew,Lmax,inv,Lrange);

    varns={dataout};
    varargout=varns(1:nargout);



elseif strcmp(coefin,'demo1')
    % Evaluate topography over North America
    [rad,theta,phi]=randinpolyradius('namerica',400,(6371+400)*1000,100000);
    lon = phi*180/pi;
    lat = 90-theta*180/pi;

    thetopo = fralmanac('EGM2008_Topography','SHM');
    Lmax=180;
    thetopo = thetopo(1:length(addmon(Lmax)),:);

    dataout = potupxyz(thetopo,6371*1000,lon,lat,rad,Lmax);

    XY = namerica();
    figure
    plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
    hold on
    scatter(lon,lat,10,dataout,'filled')

end


