function diferm(a,b,tolex,sev)
% DIFERM(a,b,tolex,sev)
% 
% The mute form of DIFER, i.e. difer(a-b,tolex,sev,NaN)
%
% a         The first vector offered for comparison
% b         The second vector offered for comparison
% tolex     Tolerance exponent (default: 10 for 1e-10)
% sev       0 produces WARNING upon failure [default]
%           1 produces ERROR upon failure 
%           2 invokes KEYBOARD upon failure
%
% Last modified by fjsimons-at-alum.mit.edu, 03/20/2020

% Any second input
defval('b',0)
% The tolerance exponent
defval('tolex',[])
% The severity index
defval('sev',[])

% Run the test
difer(a-b,tolex,sev,NaN)
