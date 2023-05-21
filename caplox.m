function [lon2,lat2]=caplox(lonlat,TH,N,method)
% [lon2,lat2]=CAPLOX([lon1 lat1],TH,N)
%
% Calculates the location of a circle on the sphere crossed with the
% outlines of the continents.
%
% INPUT:
%
% [lon1 lat1]      Longitude, latitude of center(s) [degrees]
% TH               Radius, in degrees [default: 15]
% N                Number of points [default: 100]
%
% OUPUT:
%
% [lon2,lat2]      The requested points [degrees or Mollweide]
%
% See also: LATITUDE, LONGITUDE, CAPLOC
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2021

defval('lonlat',[-25 10])
defval('TH',65)

defval('N',100)
defval('method',1)

% Calculate the circle
angl=linspace(0,360,N);
delta=TH*fralmanac('DegDis','Earth')/1000;
[lon2,lat2]=grcazim(lonlat,delta,angl);
lon2(lon2>180)=lon2(lon2>180)-360;

% Centered on the meridian?
conz=1;
if conz==0
  % Get the continents - like in PLOTCONT
  defval('ddir',fullfile(getenv('IFILES'),'COASTS'))
  fid=fopen(fullfile(ddir,'cont.mtl'),'r','b');
  cont=fread(fid,[5217 2],'uint16');
  fclose(fid);
  cont=cont/100-90;
  cont(cont==max(max(cont)))=NaN;
else
  [~,cont]=maprotate([],[]);
end

% Begin by finding all continental points that are inside the circle
in=inpolygon(cont(:,1),cont(:,2),lon2,lat2);
cont(~in,:)=NaN;

% Then for each azimuth find the points and take the inner radius?
% Or calculate all azimuths for all points, bin them

% Edit out England? Spuriously removes other things
%XY=england; XY(XY(:,1)>180,1)=XY(XY(:,1)>180,1)-360;
%[C,Ic,IX]=intersect(cont,XY,'rows');

% Remove, e.g., Iceland by excluding handdrawn boxes?
% This works for TH=65 but has to be handedited
ax=[-24.4289  -13.2650   61.2055   68.6062;...
    -84.6932  -67.6656   16.3281   23.7487;...
    -12.0142    1.9796   50.7961   60.0728;...
    -5.9774    1.6643    49.9282   52.1100;...
    -59.3144  -52.1811   46.1328   49.8000;...
    -55.8137  -55.4431   51.1717   51.4883;...
    -57.5438  -55.2695   49.2515   51.1949;...
   -104.6383  -69.0724  -40.5857    6.8965;...
    -91.1963  -84.0030    8.5333   13.7169;...
    -91.4466  -89.3057   13.2858   14.8286;...
    -0.4969   42.4864   13.8457   44.8200;...
     8.9328   23.4210   38.7351   49.1754;...
    -5.2156   -0.4692   34.3656   37.7859;...
     9.2356   25.8408   48.8525   60.8184;...
    17.0586   30.5769   54.0010   70.7644;...
   -79.9360  -78.1397    6.8998    9.1272;...
   -84.0833  -82.8935    8.3059    9.7813;...
    -82.8963  -80.0105    5.9117   8.4902;...
    -86.8801  -68.5826   49.1935   71.8834];

for index=1:size(ax,1)
  ai=[cont(:,1)>ax(index,1) & cont(:,1)<ax(index,2)] & ...
     [cont(:,2)>ax(index,3) & cont(:,2)<ax(index,4)];
  cont(ai,:)=NaN;
end



figure(gcf)
clf
pc=plot(cont(:,1),cont(:,2));
hold on
ps=plot(lon2,lat2);
hold off
axis equal


% hold on
% % Get the continents in viewed from the center of the circle
% [a,b,XYZ]=plotcont([],[],11,[],[],lonlat);
% keyboard
% p=plot(lonlat(1),lonlat(2))
% hold off

