function deggies(ah,varargin)
% DEGGIES(ah)
% DEGGIES(ah,1) % Only x
% DEGGIES(ah,2) % Only y
%
% This function puts degree signs on the ticklabels of the
% axes with handle 'ah'. Takes vectors!
%
% See also DEGS
%
% Last modified by fjsimons-at-alum.mit.edu, 12/26/2006

defval('ah',gca)

if nargin>1
  swit=varargin{1};
else 
  swit=3;
end

for index=1:length(ah)
  XTL=get(ah(index),'XTickLabel');
  YTL=get(ah(index),'YTickLabel');

  space=32;
  circ=176;

  % This is for new Matlab, which renders badly in postscript
  XTL=[XTL repmat(str2mat(space),size(XTL,1),1)];
  YTL=[YTL repmat(str2mat(space),size(YTL,1),1)];

  mnx=size(XTL);
  mny=size(YTL);

  indx=sum(XTL==' ',2);
  indy=sum(YTL==' ',2);

  extx=(circ-space)*(sparse(1:mnx(1),mnx(2)-indx+1,1));
  exty=(circ-space)*(sparse(1:mny(1),mny(2)-indy+1,1));

  if ~isempty(XTL) && swit~=2
    XTL=str2mat(abs(XTL)+extx); 
    % For latest Matlab; after this, do $UFILES/degs
    % Actually, this doesn't help
    % XTL=[repmat(' ',size(XTL,1),1) XTL];
    % HOW ABOUT setstr(hex2dec('B0')) THIS DOESN'T HELP EITHER
  end
  if ~isempty(YTL) && swit~=1 
    YTL=str2mat(abs(YTL)+exty); 
    % For latest Matlab; after this, do $UFILES/degs
    % Actually, this doesn't help
    % YTL=[repmat(' ',size(YTL,1),1) YTL];
  end
  
  switch swit
    case 1
      set(ah(index),'XTickLabel',XTL)
    case 2
      set(ah(index),'YTickLabel',YTL)
    otherwise
      set(ah(index),'XTickLabel',XTL)
      set(ah(index),'YTickLabel',YTL)
  end
end
