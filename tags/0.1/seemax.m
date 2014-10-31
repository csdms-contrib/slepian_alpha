function seemax(ah,wat)
% function SEEMAX(ah,xycz)
%
% Equalizes axis ranges between plots. 
%
% INPUT: 
% 
% ah        Axis handles
% xycz      One or many of:
%           1 Equalizes xlim
%           2 Equalizes ylim
%           3 Equalizes clim
%           4 Equalizes zlim
%
% EXAMPLE: 
%
% clf
% a=subplot(222); plot(randn(10,3))
% b=subplot(221); plot(randn(10,6)*5)
% seemax([a b],2)
%
% Last modified by fjsimons-at-alum.mit.edu, January 15th, 2003

prop1={'Xlim','Ylim','Clim','Zlim'};

% Capitals in these positions matter
prop2={'XData','YData','','ZData'};

for index=1:length(wat)
  wot=wat(index);
  counter=0;
  switch wot 
   case 3
    for indo=1:length(ah)
      counter=counter+1;
      mindata(counter)=min(get(ah(indo),prop1{wot}));
      maxdata(counter)=max(get(ah(indo),prop1{wot}));
    end
   otherwise
    for indo1=1:length(ah)
      kids=get(ah(indo1),'Children');
      for indo2=1:length(kids)
	counter=counter+1;
	if isfield(get(kids(indo2)),prop2{wot})
	  mindata(counter)=min(min(get(kids(indo2),prop2{wot})));
	  maxdata(counter)=max(max(get(kids(indo2),prop2{wot})));
	end
      end 
    end
  end
  for indo1=1:length(ah)
    set(ah(indo1),prop1{wot},[min(mindata) max(maxdata)])
  end
end

