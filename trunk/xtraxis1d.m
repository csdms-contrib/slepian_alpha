function varargout=xtraxis1d(aha)
% XTRAXIS1D(aha)
% xtra=XTRAXIS1D(aha)
%
% Creates wavelength axis when the original axis was in radians per km.
% For one-dimensional logarithmic plots with positive frequencies only.
%
% See also XTRAXIS1DLIN, XTRAXIS2D, XTRAXIS
%
% Watch out: need to manually adjust xticks sometimes
%
% Last modified by fjsimons-at-alum.mit.edu, October 21st 2003

xel=xlim(aha);
yel=ylim(aha); yti=get(aha,'Ytick');
xtra=axes('Position',get(aha,'Position'));
set(aha,'box','off')
set(xtra,'XAxisLocation','top','YAxisLocation','right')
set(xtra,'Ylim',yel,'YTick',yti,'YTickLabel',[])
set(xtra,'Xlim',sort(2*pi./xel),'XTick',[100:100:1000 2000])
set(xtra,'Color','none','XDir','reverse','FontS',get(aha,'FontS'),...
    'XTickLabel',{'100' ' ' '300' ' ' ' ' '600' ' ' ' ' ' ' '1000' '2000'})
set(xtra,'XScale','log','xminortick','off')

if nargout
  varargout{1}=xtra;
end

