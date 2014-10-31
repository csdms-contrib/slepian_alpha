function cols=halverange(data,fax,fidi)
% cols=HALVERANGE(data,fax,fidi)
%
% Makes a SYMMETRIC half-range color limit... or set it to some
% percentage as referred to max(abs(data))
%
% INPUT:
%
% data       The data vector or matrix
% fax        50, for 50% [default]
% fidi       NaN do not display the message
%
% Last modified by fjsimons-at-alum.mit.edu, 3/15/2011

defval('fax',50)
defval('fidi',1)

cols=[-max(abs(data(:))) max(abs(data(:)))]*fax/100;

if ~isnan(fidi)
  disp(sprintf('Range is %8.4e to %8.4e',cols(1),cols(2)))
end

