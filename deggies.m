function deggies(ah,swit)
% DEGGIES(ah,swit)
%
% Puts degree signs on axis ticklabels
%
% INPUT:
%
% ah       The axis handle(s)
% swit     1 Only on the x-axis
%          2 Only on the y-axis
%          3 On both the x and the y axes (default)
%
% See also DEGS, a Unix script to fix a PostScript error
% 
% Tested on 8.3.0.532 (R2014a)
%
% Last modified by fjsimons-at-alum.mit.edu, 03/12/2020

defval('ah',gca)
defval('swit',3)

% Apply to all the axis handles necessary
for index=1:length(ah)
  XTL=get(ah(index),'XTickLabel');
  YTL=get(ah(index),'YTickLabel');

  % ASCII codes for special things
  space=32;
  circ=176;

  if verLessThan('matlab', '9.0.0')
    % This is what renders badly in PostScript
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
      % Actually, this doesn't help
      % XTL=[repmat(' ',size(XTL,1),1) XTL];
      % HOW ABOUT setstr(hex2dec('B0')) THIS DOESN'T HELP EITHER
    end
    if ~isempty(YTL) && swit~=1 
      YTL=str2mat(abs(YTL)+exty); 
      % Actually, this doesn't help
      % YTL=[repmat(' ',size(YTL,1),1) YTL];
    end
  else
    if ~isempty(XTL) && swit~=2
      for ondex=1:size(XTL,1)
        XTL{ondex}=sprintf('%s%s',XTL{ondex},circ);
      end
    end
    if ~isempty(YTL) && swit~=1 
      for ondex=1:size(YTL,1)
        YTL{ondex}=sprintf('%s%s',YTL{ondex},circ);
      end
    end
  end
  % Effectuate the switch
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
