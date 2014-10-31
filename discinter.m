function [newval,req,dispos,avma]=discinter(pos,val,req)
% [newval,req,dispos,avma]=DISCINTER(pos,val,req)
%
% For non-monotonous vectors (such as with discontinuities, specified 
% by two or three times the same value), interpolates, i.e. finds the
% interpolated value of input values.
% Returns different values than requested!
%
% Returns NaN outside the range, as in regular interpolation.
% Also for regular interpolation, i.e. without discontinuities.
% 'pos' must not start with a discontinuity
%
% The 'newval' belong to an upgoing sorted 'req', possibly with 
% fewer entries than the input 'req', if the latter was discontinuous
% either as 2p or 3p.
%
% Discontinuities at zero are treated as normal ones, so they might 
% end up having negative values. So when doing inner products
% and stuff, make sure to specify values a little bit out of the range.
%
% Everything discontinuity related will get differences (intervals)
% of 0.5 and 1.5 times smint10, all the rest are real intervals.
% Changed smint/10 to smallest of 1 m or smint/10.
%
% Also returns the vector 'dispos' from DISCIDENT(POS)
% Also returns the averaging matrix 'avma' from DISCIDENT(POS)
%
% See DISCINTEREX, and NODINTER
% 
% For example, an request vector [15 30 80 140 200 300 400 400] has
% a discontinuity at 400. But the position/value pair has discons at
% 15 80 and 400, specified as a repeated value (2p).
% The one at 400 is common, but 15 and 80 aren't
% The result is: you get values at
% [14 14.5 16 30 79 80.5 81 140 200 300 399 401]
% If you ask what the value at 15 is, you can't get a straight answer, but
% get three values instead. 
% If you ask what the value is at 400 and 400, you get the values at 399 and 401
%
% Last updated by fjsimons-at-alum.mit.edu, October 18th, 2002

[pos,val,req]=deal(pos(:),val(:),req(:));
if pos(1)>pos(end)
  pos=flipud(pos);
  val=flipud(val);
  fli=1;
else
  fli=0;
end

% If requested values are 2p-3p-discontinuities in the request array, then
% destroy, not even leaving the middle value intact for 3p.
difreq=diff(req);
smint=min(difreq(~~difreq));
smint10=min(1,smint/10);
if sum(~difreq)
  [a,b]=discident(req);
  % Leave quasi-3p for quadrature; work with NaN's
  % Smallest interval of req
  req(~isnan(a))=req(~isnan(a))-(smint10)*a(~isnan(a));
  % For 3p, adjust the middle value, for 2p, leave intact
  req(~b)=req(~b)-(smint10)/2;
end

% This makes [400 400] into [399 401]

req=sort(req);

outsid=req<min(pos) | req>max(pos);
reqout=req(outsid);
req=req(~outsid);
outvol=repmat(NaN,sum(outsid),1);

% What are the discontinuities of 'pos'?
[c,d,posdisc,avma]=discident(pos);
dispos=c;
% If it was flipped, unflip for output
if fli==1
  dispos=flipud(dispos);
  avma=flipud(avma);
end

% So now the requested vector 'req' is free of discontinuities.
% Now need to make sure that the values of 'req'
% that would be equal to an element of 'pos' are treated
% properly. But 'req' can still have values equal to a discontinuity in pos
% this was the original point of this program.
flag1=0;
if ~isempty(intersect(pos,req))
  [pdi,inp,inr]=intersect(posdisc,req);
  if ~isempty(pdi)
    % Remediate in the request vector, make it want a "3p-discontinuity"
    % 3p for integration purposes
    req(inr)=pdi(:)-repmat(smint10,length(pdi),1);
    req=[req ; pdi(:)+repmat(smint10,length(pdi),1)];
    req=[req ; pdi(:)+repmat(smint10/2,length(pdi),1)];
    req=sort(req);
  end
  flag1=1;
  % Else you treat those values
  % I can take them out, knowing they are not discs in either
  % vector, lookup their unambiguous value and put them in at 
  % the end, knowing the output was sorted anyway.
  [reqvol,posinpos,posinreq]=intersect(pos,req);
  newvol=val(posinpos);
  req=req(~ismember(1:length(req),posinreq));
end

% If there are 2p or 3p discontinuities in the positions
if sum(~diff(pos))
  posorg=pos;
  % This removes the discontinuities from 'pos';
  % it is assumed that 'pos' and 'req' have no common members.
  % This is where 'req' gets sorted!
  pos=unique([pos ; req]);
  ins=[1:length(pos)]';
  
  % reqpos is the index of req in [pos-disc+req]
  reqpos=ins(ismember(pos,req));
  % reqpos is index of one above req in [pos-disc]
  reqpos=reqpos-[0:length(reqpos)-1]';

  % This need to be the position shift needed to translate
  % the original pos index into pos-disc 
  [a,ind,c]=degamini(posorg);
  ind=cumsum([0 ind-1]); ind=ind(:);
  
  % What is the index, in the original value matrix, of the
  % parameter just down the one we are requesting?
  % And then ups-1 is the one right below
  ups=reqpos+ind(reqpos);
  newval=val(ups-1)+(req-posorg(ups-1))./...
      (posorg(ups)-posorg(ups-1)).*(val(ups)-val(ups-1));
else
  newval=interp1(pos,val,req,'linear');
end

if flag1==1
  outrows=sortrows([req newval ; reqvol newvol ; reqout  outvol]);
else
    outrows=sortrows([req newval ; reqout  outvol]);
end
[req,newval]=deal(outrows(:,1),outrows(:,2));
