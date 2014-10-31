function [g,s,h,hh]=row2stats(col1,col2,ybine)
% [g,s,h,hh]=ROW2STATS(col1,col2,ybine)
%
% Find statistics of data presented belonging to different groups
%
% INPUT:
%
% col1    The first column of data of which groups are made 
% col2    The second column of data 
% ybine   The bin edges to make histograms of the second data column
%
% OUTPUT:
%
% g       The different groups identified in the first data column
% s       Structure array with mean, median, etc of data per group 
% h       Matrix with histogram counts per group in the bins provided
% hh      Matrix with histogram counts globally normalized
% 
% See also BINDENS 
%
% Last modified by fjsimons-at-alum.mit.edu, 12/21/2011

defval('nybins',10)
defval('ybin',range(col2)/nybins);
defval('ybine',linspace(min(col2),max(col2),nybins+1));

% Sort the data in the first column
[col1,I]=sort(col1,1);
% And have the second column follow
col2=col2(I);
% The beginning index of every new value in column 1
beg=[1 find(diff([col1(1) col1(:)']))];
% The ending index of every separate magnitude
ent=[beg(2:end)-1 length(col1)];

% Now collect statistics of the values in column 2
for index=1:length(beg)
  coli=beg(index):ent(index);
  % Identify the category in the first column
  g(index)=col1(beg(index));
  % Identify the element of the second column
  colrel=col2(coli);
  % Construct the statistics of the second column
  s.mean(index)=nanmean(colrel);
  s.median(index)=nanmedian(colrel);
  s.nonnans(index)=sum(~isnan(colrel));
  s.variance(index)=nanvar(colrel);
  s.p25(index)=prctile(colrel,25);
  s.p75(index)=prctile(colrel,75);
  if nargout>2
    h(:,index)=histc(colrel,ybine);
    hh(:,index)=h(:,index);
    h(:,index)=h(:,index)/sum(h(:,index));
  end
end

if nargout>2
  % Now remember the oddity about the HISTC bins
  h(end-1,:)=h(end-1,:)+h(end,:);
  h=h(1:end-1,:);
  hh(end-1,:)=hh(end-1,:)+hh(end,:);
  hh=hh(1:end-1,:);
end

