function varargout=krijetem(pols,pos)
% [ah,ha,H]=KRIJETEM(pols,pos)
%
% Creates subplots returns handles to them.
%
% INPUT:
%
% pols    Panel specifiers, e.g. [221 224], or ['2,2,1';'2,2,4']
%         OR: lbwh position coordinates in case the next variable is 'pos'
% pos     Nothing [default] or 'pos' if the first input contains coordinates
%
% OUTPUT:
%
% ah      The handles in horizontal running order
% ha      The handles in vertical running order (only when panel numbers
%         are complete sequences as provided by, e.g. SUBNUM)
% H       The handles in matrix form (once again only when this makes
%         sense), which is useful for passing them onto SERRE
%
% SEE ALSO: SUBNUM, SERRE
%
% [ah,ha,H]=krijetem(subnum(3,3));
% [ah,ha]=krijetem(['2,2,3' ; '3,3,4'])
%
% Last modified by fjsimons-at-alum.mit.edu, 04/03/2009

% Supply default
defval('pos',[])

if ~isempty(pos) & strcmp(pos,'pos')
  % Specify actual positions
  for index=1:size(pols,1)
    ah(index)=axes('position',pols(index,:));    
  end
else
  % Specify subplot panel indices as concatenated numbers - only through 339
  if ~isstr(pols)
    pols=pols(:)';
    % It's three numbers and each of them is distinct
    nwy=floor(pols(1)/100);
    nwx=floor((pols(1)-floor(pols(1)/100)*100)/10);
    for index=1:length(pols)
      ah(index)=subplot(pols(index));
    end
  else
    % Specify subplot panel indices as a string of comma-separated numbers
    [nwy,cpl]=strtok(pols(1,:),','); 
    nwy=eval(nwy);
    nwx=eval(strtok(cpl(1,:),','));
    for index=1:size(pols,1)
      eval([ 'ah(index)=subplot(',pols(index,:),');'])
    end
  end
  if length(pols)==nwx*nwy
    % Only here does the column/row matrix ordering make any sense
    H=reshape(1:nwy*nwx,nwx,nwy)';
    ha=ah(H(:));
    H=ah(H);
  end
end

% If you haven't got them by now, don't bother
defval('ha',ah)
defval('H',ah)

% Provide desired output
varns={ah,ha,H};
varargout=varns(1:nargout);
