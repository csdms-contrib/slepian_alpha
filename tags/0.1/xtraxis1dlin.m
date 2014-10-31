function xtra=xtraxis1dlin(aha,elst,elstm,faxt)
% xtra=XTRAXIS1DLIN(aha,elst,elstm,faxt)
%
% Creates wavelength axis when the original axis was in radians per km.
% For one-dimensional linear plots with positive frequencies only.
%
% INPUT:
%
% aha       Axis handle (default: gca)
% elst      Which values on the new axis to annotate
% elstm     And a cell array with their labels in reverse
% faxt      Factor multiplying the frequency (default: 2*pi)
%
% OUTPUT:
%
% xtra      Handle to the new axis
%
% See also XTRAXIS, XTRAXIS1D, XTRAXIS2D
%
% Last modified by fjsimons-at-alum.mit.edu, 06.10.2005

% Ready for which ticks and ticklabels?
defval('elst',[100:100:1000 2000]);
defval('elstm',{ '2000' ' ' ' ' ' ' '700' ' ' ' ' '400' '300' '200' '100'});
defval('faxt',2*pi);

% Which ones to annotate?
annot=sort(faxt./elst);

% Get the limits of the axis in question
xel=xlim(aha);
yel=ylim(aha); 
yti=get(aha,'Ytick');

% Create new axis on old position; not with laxis
xtra=axes('Position',getpos(aha));

% Put them on
set(xtra,'Color','none','Ylim',yel,'YTick',yti,'YTickLabel',[],...
    'Xlim',xel,'XaxisL','top','YaxisL','right','Color','none',...
    'FontSize',get(aha,'FontSize'),...
    'XTick',annot,'XTickL',elstm)

set(aha,'box','off')
