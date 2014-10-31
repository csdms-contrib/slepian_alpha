function [out,ind]=skip(in,pos)
% [out,ind]=SKIP(in,pos)
%
% Returns a vector with certain positions skipped.
%
% INPUT:
%
% in       The input vector
% pos      The indices to be skipped
%
% OUTPUT:
%
% out      The output vector
% ind      Logical array returning the skipped positions of the input
%
% Last modified by fjsimons-at-alum.mit.edu, 03/21/2012

in=in(:)';
pos=pos(:)';
pos=unique(pos);

ind=ones(length(in),1); 
ind(pos)=0; 
out=in(~~ind);
