function varargout=movev(fax,amt)
% [new,old]=MOVEV(fax,amt)
%
% Moves an axis or object handle VERTICALLY by some amount.
% 
% INPUT:
%
% fax      Axis handles or non-axis object handles
% amt      Amount of move (+ up, - down) [default: 0.1]
%
% OUTPUT:
%
% new      The new positions
% old      The old positions
%
% SEE ALSO:
%
% MOVEH
%
% Last modified by fjsimons-at-alum.mit.edu, 04/02/2009

defval('amt',0.1)

% Figure out the type of object
tip=strcmp('patch',get(fax,'type'));

% Make into a vector and do for all
fax=fax(:);

% Do this for all of the axis handles
for index=1:length(fax)
  switch all(tip)
    case 0
     % Record the current (old) position
     oldpos=getpos(fax(index));
     % Set to the new position
     newpos=repmat(0,1,length(oldpos));
     newpos(2)=amt;
     newpos=oldpos+newpos;
     set(fax(index),'Position',newpos)
   case 1
    % Record the current (old) position
     oldpos=get(fax(index),'Vertices');
     % Set to the new position
     newpos=repmat(0,size(oldpos));
     newpos(:,2)=amt;
     newpos=oldpos+newpos;
     set(fax(index),'Vertices',newpos)
  end
end

% Record the current (new) position
switch all(tip)
 case 0
  newpos=getpos(fax);
 case 1
  newpos=get(fax,'vertices');
end

% Provide desired output
varns={newpos,oldpos};
varargout=varns(1:nargout);
  
