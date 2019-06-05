function stronk=parse(strink,sepor)
% stronk=PARSE(strink,sepor)
%
% Makes a string matrix from a long string delimited by a certain
% separator and strips it of rows of blanks as well
%
% INPUT:
%
% strink     The string that you want parsed
% sepor      The single separator to guide the parsing [default: a newline]
%
% EXAMPLE:
%
% files=parse(ls(['x200' '*']));
%
% Last modified by fjsimons-at-alum.mit.edu, 06/05/2019

% Default separator
defval('sepor',sprintf('\n'))

% Go find the separator!
ent=findstr(strink,sepor)-1;

if ~isempty(ent)
  % If the string does not END on the delimiter...
  if ent(end)~=[length(strink)-1]
    ent=[ent length(strink)];
  end
  % Where the parse segments begin
  beg=[1 ent+2];

  % Initialize
  stronk= ' ';
  % Iterate
  for index=1:length(beg)-1
    % Remember STR2MAT was the predecessor to CHAR
    try
      stronk=char(stronk,strink(beg(index):ent(index)));
    catch
      stronk=str2mat(stronk,strink(beg(index):ent(index)));
    end
  end
  stronk=stronk(2:end,:);
else
  %  ... otherwise just take the whole thing
  stronk=strink;
end

% Remove rows of blanks
delt=1:size(stronk,1);
for index=1:size(stronk,1)
  if all(abs(stronk(index,:))==32)
    delt(index)=0;
  end
end

% The end result!
stronk=stronk(~~delt,:);

