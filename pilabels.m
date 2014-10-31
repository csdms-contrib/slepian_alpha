function pilabels(ah)
% PILABELS(ah)
%
% Makes X-axis labels go from 0 to pi
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 12th, 2004

for index=1:length(ah)
  xti=linspace(0,pi,5);
  xtil={'0' 'p/4' 'p/2' '3p/4' 'p'};
  
  set(ah(index),'xtick',xti,'xtickl',xtil) 
  set(ah(index),'FontN','symbol')
end


