function RGB=getcol(symb)
% RGB=GETCOL(symb)
%
% Returns RGB values corresponding to color symbols
%
% INPUT:
%
% 'symb'   A column or row vector or a cell array of single-character
%          color strings from the list y m c r g b w k
%
% OUTPUT:
% 
% RGB      A three-column matrix with RGB values
% 
% EXAMPLE:
%
% getcol([ 'y' 'm' 'k' 'g' 'y'])
% getcol([ 'y' ,'m' ,'k' ,'g'])
% getcol([ 'ymkg'])
% getcol([{'y'} ,{'m'} ,{'k'} ,{'g'}]))
%
%
% Last modified by fjsimons-at-mit.edu, 05/26/2021

% Input parsing
if ~iscell(symb)
  symb=symb(:)';
else
  symb=cat(1,symb{:})';
end

% Character mapping
symb=abs(symb)-97;

% Color identification from an alphabet from b to y
coleq=ones(24,3);
coleq(1,:) =[0 0 1]; %b
coleq(2,:) =[0 1 1]; %c
coleq(6,:) =[0 1 0]; %g
coleq(10,:)=[0 0 0]; %k
coleq(12,:)=[1 0 1]; %m
coleq(17,:)=[1 0 0]; %r
coleq(22,:)=[1 1 1]; %w
coleq(24,:)=[1 1 0]; %y 

% Color selection assignment
RGB=coleq(symb,:);
