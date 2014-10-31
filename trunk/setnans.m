function [data,thre]=setnans(data,fax)
% [data,thre]=SETNANS(data,fax)
%
% INPUT:
%
% data     The data matrix
% fax      faxth fraction of max(abs(data)) below which data are replaced
%          by NaN [default: 100]
% 
% OUTPUT:
%
% data     The new data matrix
% thre     The threshold used 
%
% Last modified by fjsimons-at-alum.mit.edu, 12/18/2012

defval('fax',100)
% Check for -Inf or +Inf and get rid of those 
isi=isinf(data);
if sum(isi(:))>0
  fprintf(1,'%i infinite values replaced by 0 and thus by NaN\n',...
          sum(isi(:)))
  data(isinf(data))=0;
end
% Threshold
thre=max(abs(data(:)))/fax;
data(abs(data)<thre)=NaN;
