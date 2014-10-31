function lmcosif=plmfilt(lmcosi,L,wlen)
% lmcosif=plmfilt(lmcosi,L,wlen)
%
% A very simple filtering routine for spherical harmonics
%
% INPUT:
%
% lmcosi    Standard [degree order cosine sine] coefficients, not
%           necessarily starting at zero degree
% L         Single degree with lowpass degree of halfway point
%           or two degrees for a bandpass filter, NaN for nothing 
% wlen      The odd length of the transition band [default: 5]
%
% OUTPUT:
%
% lmcosid   Filtered equivalent of the input array
%
% EXAMPLE:
%
% lmcosi=plm2rnd(30,-2);
% lmcosif=plmfilt(lmcosi,[7 20],5);
%% Compare the input/output ratios - knowing the degree ranges need to correspond
% [lmcosif(:,1:2) lmcosif(:,3:4)./lmcosi(addmup(lmcosif(1,1)-1)+1:addmup(lmcosif(end,1)),3:4)]
%
% Thanks to rmr-at-arizona.edu for pointing out an earlier error
% Last modified by fjsimons-at-alum.mit.edu, 08/05/2014


% Specify defaults
defval('wlen',5)
defval('L',180)

if isnan(L)
  lmcosif=lmcosi;
  return
end

% Make sure it's sorted descending
L=sort(L,'descend');

if ~mod(wlen,2)
  error('Total filter length must be odd')
end

% What is the half length of the filter
whalf=(wlen-1)/2;

% Generate the window as filter coefficients
[w,wl,wr]=fhanning(2*wlen);

% Initialize
lmcosif=lmcosi;

% Where does it all start, that determines the offset
lofs=addmup(lmcosi(1)-1);

% Repeat them the right number of times for the center degree
for index=1:length(L)
  if L(index)+whalf>lmcosi(end,1) | L(index)-whalf<lmcosi(1,1)
    error('Filter too long or bandlimit too high/low for input')
  end
  degfilt=L(index)+[-whalf:whalf];
  % Cut to the right
  if index==1 
    wrr=gamini(wr,degfilt+1)';
  end
  % Cut to the left if there is a LOWER degree coming
  % Note that FLIPUD and GAMINI do not commute!
  if index==2
    wrr=gamini(wl,degfilt+1)';
  end
  disp(sprintf('Filtering degrees %i - (%i) - %i',...
	       degfilt(1),L(index),degfilt(end)))
  % Perform the filtering 
  degind=[addmup(degfilt(1)-1)+1-lofs:...
	  addmup(degfilt(end))-lofs]';
  lmcosif(degind,3)=lmcosif(degind,3).*wrr;
  lmcosif(degind,4)=lmcosif(degind,4).*wrr;
end

% Truncate to the new maximal resolution
lmcosif=lmcosif(1:addmup(L(1)+whalf)-lofs,:);

% Truncate from the new minimal resolution
if length(L)==2
  lmcosif=lmcosif(addmup(L(2)-whalf-1)-lofs+1:end,:);
end
