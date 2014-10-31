function RGB=getcol(symb)
% RGB=GETCOL(symb)
%
% Gets RGB values of M color symbols
% 'symb' can be a vector 1XM or MX1
% or a cell structure 1XM or MX1
%
% RGB is MX3
%
% symb=[ 'y' 'm' 'k' 'g' 'y']
% symb=[ 'y' ,'m' ,'k' ,'g']
% symb=[ 'ymkg']
% symb=[{'y'} ,{'m'} ,{'k'} ,{'g'}]
%
% Posisble: ymcrgbwk

% Last modified by fjsimons-at-mit.edu, June 6th 2000
if ~iscell(symb)
  symb=symb(:)';
else
  symb=cat(1,symb{:})';
end

symb=abs(symb)-97;

coleq=ones(24,3);
coleq(24,:)=[1 1 0]; %y 
coleq(12,:)=[1 0 1]; %m
coleq(2,:)=[0 1 1]; %c
coleq(17,:)=[1 0 0]; %r
coleq(6,:)=[0 1 0]; %g
coleq(1,:)=[0 0 1]; %b
coleq(22,:)=[1 1 1]; %w
coleq(10,:)=[0 0 0]; %b

RGB=coleq(symb,:);
