function lcd=even(x)
% lcd=EVEN(x)
%
% If x is a number, returns lcd for x=length
%
% Returns LOGICAL downsampling array isolating the
% even-numbered coefficients, where the
% first sample is ODD (Matlab convention)
%
% Don't confuse with Matlab's function ODD
%
% Note that the number of 
% EVEN entries is always floor(length(x)/2)
% ODD  entries is always ceil(length(x)/2)
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 14th, 2003

if length(x)>1
  lcd=~logical(1/2-1/2*(-1).^(1:length(x)))';
else
  lcd=~logical(1/2-1/2*(-1).^(1:x))';
end



