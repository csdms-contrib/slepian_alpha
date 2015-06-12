function varargout=struct2str(strux,fid)
% strux=STRUCT2STR(strux,fid)
% STRUCT2STR(strux,fid)
%
% Prints the non-empty fields of a structure to screen or to file
%
% INPUT:
%
% strux       A structure array
%
% OUTPUT::
%
% strux       The structure array after blank removal
%
% Last modified by fjsimons-at-alum.mit.edu, 10/06/2014

defval('fid',[])

fn=fieldnames(strux);

% Get rid of the empties
for index=1:length(fn)
  struc=strux.(fn{index});
  if isempty(struc)
    strux=rmfield(strux,fn{index});
  end
end

% Reassess
fn=fieldnames(strux);

% Write out if no output requested)
if nargout==0
  if isempty(fid)
    disp(strux)
  else
    for index=1:length(fn)
      fni=fn{index};
      % fprintf(fid,'%20s:%20s\n',fni,num2str(strux.(fni)));
      fprintf(fid,'%20s:%20s\n',fni,strux.(fni));
    end
  end
end

% Optional output
varns={strux};
varargout=varns(1:nargout);
