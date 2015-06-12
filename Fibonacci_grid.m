function varargout = Fibonacci_grid(N)

% Construct Fibonacci grid, see Álvaro González, "Measurement of
% Areas on a Sphere Using Fibonacci and Latitude–Longitude Lattices",
% Mathematical Geosciences (2010) 42: 49–64


% The number of points equals: 2N + 1.

defval('N',1000);
defval('WhichAngle','ComplGoldenAngle')
%Golden ratio
Phi = (1+sqrt(5))/2;

lat = asin((2*[-N:1:N])/(2*N+1))'*(180/pi);

if strcmp(WhichAngle,'ComplGoldenAngle')
    lon = mod([-N:1:N],Phi)'*(360/Phi);
elseif strcmp(WhichAngle,'GoldenAngle')
    lon = -mod([-N:1:N],Phi^2)'*(360/Phi^2);
end

lon(lon < 0)   = lon(lon < 0) + 360;
lon(lon > 360) = lon(lon > 360) - 360;

% axesm globe
% plotm(lat,lon,'.')
% plot(lon,lat,'.')

% Prepare output
vars={lon,lat};
varargout=vars(1:nargout);

end
