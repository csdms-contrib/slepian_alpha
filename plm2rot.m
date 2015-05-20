function varargout=plm2rot(lmcosi,alp,bta,gam,method,rlcp)
% [lmcosip,spec1,spec2]=PLM2ROT(lmcosi,alp,bta,gam,method,rlcp)
%
% Rotates a scalar FIELD expanded into real spherical harmonics on the
% unit sphere surface using Euler angles in an ACTIVE rotation convention
% (the inverse of DT pp 920-924). The results are spherical-harmonic
% expansion coefficients for the new field, in the original, unrotated
% coordinate system. To invert, flip signs and orders of the Euler angles.
%
% INPUT:
%
% lmcosi           Matrix with [l m Ccos Csin] in order, m>=0
% alp, bta, gam    Euler angles [degrees], POINTS are rotated over 
%                     alp (0<360) around z increasing from y to x, then 
%                     bta (0<180) around old y increasing from x to z, then
%                     gam (0<360) around old z increasing from y to x.
%                  Equivalently but alternatively, POINTS are rotated over
%                     gam (0<360) around z increasing from y to x, then
%                     bta (0<180) around new y increasing from x to z, then
%                     alp (0<360) around new z increasing from y to x.
%                  Equivalently but alternatively, POINTS are rotated
%                     gam-pi/2 around z increasing from y to x, then
%                         pi/2 around new y increasing from x to z, then
%                         0    around new z increasing from y to x, then
%                  followed by a second rotation of the POINTS,
%                     bta      around z increasing from y to x, then
%                        -pi/2 around new y increasing from x to z, then
%                     alf+pi/2 around new z increasing from y to x.
%                  Note that in order to do perform this ACTIVE rotation
%                  as described above, we stick to the PASSIVE rotation
%                  framework described by Dahlen and Tromp, App. C.8., in
%                  other words the equations for PASSIVE are the ones for
%                  ACTIVE but backwards, flipping their signs and
%                  orderings. Hence the perversity of writing about
%                  rotating points and finding the expansion in the
%                  unrotated system but talking about "new" axes in the
%                  description... it's AS IF the original axes were
%                  rotated like the "points", which we use strictly to
%                  clarify the operations...
%                  Find the appropriate rotation angles for simple
%                  geographic location, that is, the former North Pole moves:
%                  alp=0;        % Around original z axis, "clockwise"
%                  bta=lat-90;   % To desired latitude, around old y axis
%                  gam=180-lon;  % To desired longitude, around old z axis 
%                  which should be consistent with the arguments in PTOSLEP
%
% method           'dlmb' using DLMB with bta=90 [default] after decomposition
%                  'blanco' using BLANCO with bta=90 after decomposition
%                  'blanco2' using BLANCO without the rotation decomposition
% rlcp             1 if coefficients belong to REAL harmonics [default]
%                  2 if coefficients belong to COMPLEX harmonics 
%
% OUTPUT:
%
% lmcosip          Matrix with [l m Ccosp Csinp] rotated
% spec1            Spectral density of input field
% spec2            Spectral density of rotated field (identical to precision)
%
% EXAMPLES:
%
% plm2rot('demo1') % For a random input with some spectral slope
% plm2rot('demo2') % For Master's example coefficient set
% plm2rot('demo3') % For Wieczorek-Simons taper, nice figure
% plm2rot('demo4') % For a random pure spherical harmonic
% plm2rot('demo5') % For a random azimuthal cap, nice figure
% plm2rot('demo6') % For a random non-azimuthal cap, nice figure
% plm2rot('demo7') % Difference between BLANCO and DLMB
% plm2rot('demo8') % Presentation-quality figure
%
% SEE ALSO: PTOSLEP, EQPOTENTIAL, SDWTARGET
%
% Last modified by fjsimons-at-alum.mit.edu, 05/23/2013

% SHOULD BUILD IN SPECIAL CASES FOR 0 ANGLES, AS IN GEOBOXCAP

if ~isstr(lmcosi)
  defval('alp',0)
  defval('bta',0)
  defval('gam',0)
  defval('method','dlmb')
  defval('rlcp',1)
  defval('xver',0)
  
  if alp==0 && bta==0 && gam==0
    lmcosip=lmcosi;
    spec1=[]; spec2=[];
    varn={lmcosip,spec1,spec2};
    varargout=varn(1:nargout);
    return
  end
  
  % disp(sprintf('Rotation over %i %i %i to L = %i',...
  %	       round(alp),round(bta),round(gam),lmcosi(end,1)))
  
  % This amounts to the split into TWO successive rotations whereby the
  % only colatitudinal (hard) rotation is over 90 degrees (easy), twice:
  % difer(rotz(g)*roty(b)*rotz(a)-...
  %    [rotz(g+pi/2)*roty(pi/2)*rotz(0)]*[rotz(b)*roty(-pi/2)*rotz(a-pi/2)])

  % Simple conversion to radians
  alp=alp*pi/180;
  bta=bta*pi/180;
  gam=gam*pi/180;
  
  % Determine maximum L  
  L=max(lmcosi(:,1));
  % Initialize rotated array 
  lmcosip=lmcosi;

  % Methods to calculate the beta-rotation matrix
  switch method
   case 'dlmb'
    % Use Master's method, see also McEwen (2006) and elsewhere
    % disp('Using DLMB method')
    D=dlmb(L);
   case 'blanco'
    disp('Using Blanco method')
    warning('Slower, but still works')
    % Use algorithm of Blanco, Florez and Bermejo (1997).
    D=blanco(L,90);
   case 'blanco2'
    % This is a direct way... which doesn't actually work. Yet.
    % Note that I changed this by a minus sign on 3/10/2010
    [D,DC]=blanco(L,-bta*180/pi);
   otherwise
    error('Specify valid method')
  end

  if nargout>1
    spec1=plm2spec(lmcosi);
  else
    spec1=[];
  end

  if rlcp==1
    lmcosi=rsh2cpx(lmcosi);
  end
  
  % Since the algorithm is for complex harmonics

  if strmatch('blanco2',method)
    error('Not supported')
    % This is a direct method without splitting the rotation 
    % Well summarized in McEwen's 2006 paper
    for l=0:L
      [C,b,e]=shcos(lmcosi,l);
      S=shsin(lmcosi,l);
      cpl=[C ; S(2:end)];
      cplp=DC{l+1}*cpl;
      lmcosip(b:e,3)=cplp(1:l+1);
      lmcosip(b:e,4)=[0 ; cplp(l+2:end)];
    end
    if rlcp==1
      lmcosip=cpx2rsh(lmcosip);
    end
    varargout{1}=lmcosip;
    if xver | nargout>2
      spec2=plm2spec(lmcosip);
      difer(abs(sum(spec1-spec2))/L);
    else
      spec2=[];
    end
    return
  end

  % If not 'blanco2', use the decomposition and make the polar rotation
  % matrix by either blanco or dlmb for polar angle of 90 or -90.

  % warning('Changes. Updates consistent with ROTTP')
  
  % The STEPS formalism is in line with the conventions of ROTTP
  % STEP 1: PASSIVE azimuthal rotation over ALPHA-PI/2
  [C,S]=rotcof(lmcosi(:,3),lmcosi(:,4),alp-pi/2);

  % Loop over all degrees l
  cbet=cos([0:L]*bta)';
  sbet=sin([0:L]*bta)';
  for l=0:L
    li=l+1;
    % Loop over all orders m
    % Make a custom alternating matrix with zeros on the diagonal
    [i,j]=meshgrid(1:li,1:li);
    IC=mod((i+j)+mod(li,2),2)*2;  IC(:,1)=1;
    IS=~mod((i+j)+mod(li,2),2)*2; IS(:,1)=1;

    % STEP 2: PASSIVE colatitudinal rotation over -PI/2
    Cp=D{li}'.*IC*shcos(C,l);
    Sp=D{li}'.*IS*shsin(S,l);

    % STEP 3: PASSIVE azimuthal rotation over BETA
    % sign flipped as in ROTCOF
    Cpp=Cp.*cbet(1:l+1)+Sp.*sbet(1:l+1);
    Spp=Sp.*cbet(1:l+1)-Cp.*sbet(1:l+1); 

    % STEP 4: PASSIVE azimuthal rotation over ZERO
    
    % STEP 5: PASSIVE colatitudinal rotation over PI/2
    Cpp=D{li}.*IC*Cpp;
    Spp=D{li}.*IS*Spp;

    % Fill up the matrix
    b=addmup(l-1)+1;
    e=addmup(l);

    Crot(b:e,1)=Cpp;
    Srot(b:e,1)=Spp;
  end

  % STEP 6: PASSIVE azimuthal rotation over GAMMA+PI/2
  [Cp,Sp]=rotcof(Crot,Srot,gam+pi/2);
  lmcosip(:,3)=Cp;
  lmcosip(:,4)=Sp;
  if rlcp==1
    lmcosip=cpx2rsh(lmcosip);
  end
  
  if xver | nargout>2
    % Verify that spectrum remains unchanged
     spec2=plm2spec(lmcosip);
     difer(abs(sum(spec1-spec2))/L);
  else
    spec2=[];
  end

  % Provide output as needed
  varn={lmcosip,spec1,spec2};
  varargout=varn(1:nargout);
else
  demonr=lmcosi;
  switch demonr
   case 'demo1'
    lmax=ceil(20);
    [m,l,mzero]=addmon(lmax);               
    c=randn(addmup(lmax),2).*([l l].^(-1)); 
    c(mzero,2)=0; c(1,1)=1;
    lmcosi=[l m c]; 
    th0=NaN;
   case 'demo2'
    lmcm=load(fullfile(getenv('IFILES'),'masters')); 
    lmcosi=lmcm.mas; l=6; % These are complex
    lmcosi=cpx2rsh(lmcosi);
    th0=NaN;
   case 'demo3'
    th0=40;
    Lmax=36;
    [E,V,N,th,C]=sdwcap(th0,Lmax);
    [m,l,mzero]=addmon(Lmax);
    lmcosi=[l m zeros(length(l),2)];
    lmcosi(mzero,3)=C(:,1);
    th0=NaN;
   case 'demo4'
    Lmax=30;
    [em,el,mz]=addmon(Lmax); M=round(rand*Lmax);
    c=zeros(length(el),2);
    c(mz(end)+M(1),1)=1;
    lmcosi=[el em c];
    th0=NaN;
   case 'demo5'
    th0=40;
    Lmax=36;
    [E,V,N,th,C]=sdwcap(th0,Lmax);
    [m,l,mzero]=addmon(Lmax);
    lmcosi=[l m zeros(length(l),2)]; 
    lmcosi(mzero,3)=C(:,1); 
    lmcosi(mzero,4)=0;
   case 'demo6'
    th0=40;
    Lmax=36;
    [E,V,N,th,C]=sdwcap(th0,Lmax);
    [m,l,mzero]=addmon(Lmax);
    lmcosi=[l m zeros(length(l),2)];
    lmcosi(mzero+1,4)=C(:,2);
    lmcosi(mzero,4)=0;
   case 'demo7'
    lmax=ceil(20);
    [m,l,mzero]=addmon(lmax);               
    c=randn(addmup(lmax),2).*([l l].^(-1)); 
    c(mzero,2)=0; c(1,1)=1;
    lmcosi=[l m c]; 
    alp=round(rand*360);
    bta=round(rand*180);
    gam=round(rand*360);
    lmcosip1=plm2rot(lmcosi,alp,bta,gam,'dlmb');
    lmcosip2=plm2rot(lmcosi,alp,bta,gam,'blanco');
    dif=abs(sum(sum(lmcosip1(:,3:4)-lmcosip2(:,3:4))));
    disp(sprintf('BLANCO vs DLMB %8.3e',dif))
    if dif>1e-10
      error('Something wrong - correct beta angle?')
    end
    disp('No screen output')
    return
   case 'demo8'
    clf
    th0=40;
    Lmax=36;
    [E,V,N,th,C]=sdwcap(th0,Lmax);
    [m,l,mzero]=addmon(Lmax);
    alp=30;
    bta=-110;
    gam=-270;
    ah=krijetem(subnum(2,2));
    lmcosi=[l m zeros(length(l),2)];
    [xgb,ygb]=mollweide(-gam*pi/180,(bta+90)*pi/180,pi);
    [xa,ya]=longitude(-gam,3);
    [xo,yo]=latitude(bta+90,3);
    [xi,yi]=caploc([-gam bta+90],th0,[],2);
    for index=1:4
      lmcosi(mzero,3)=C(:,index);
      lmcosip=plm2rot(lmcosi,alp,bta,gam,'dlmb');
      th0=NaN;
      axes(ah(index))
      [ch(index),ph(index)]=plott3(lmcosip);
      hold on
      ps(index)=plot(xgb,ygb,'s');
      ha(index)=plot(xa,ya,'r');
      ho(index)=plot(xo,yo,'r');
      hi(index)=plot(xi,yi,'r');            
    end   
    set(ah,'CameraV',6)
    movev(ah(1:2),-.15)
    return
   case 'demo9'
    % lmcosi=rindeks(fralmanac('GTM3AR','SHM'),1:addmup(72));
    % load(fullfile(getenv('IFILES'),'TOPOGRAPHY','EARTH','GTM3AR-72'))
    lmcosi=rindeks(fralmanac('GTM3AR','SHM'),1:addmup(36));
    load(fullfile(getenv('IFILES'),'TOPOGRAPHY','EARTH','GTM3AR-36'))
    th0=NaN;
   otherwise
    error('Not a valid demo.')
  end
  clf
  % Now perform (random) rotation
  alp=round(rand*360);
  % North pole moves over bta from x to z; North Pole comes down toalp
  % center if bta=90; Antarctica comes up to center if bta=-90
  bta=round(rand*-180);
  gam=round(rand*-360);
  lmcosip1=plm2rot(lmcosi,alp,bta,gam,'dlmb');

  % And plot the stuff up
  ah=krijetem(subnum(4,2));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(1))
  [data,labs(1,:)]=plott1(lmcosi);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(2))
  [data2,labs(2,:)]=plott1(lmcosip1);
  seemax(ah(1:2),[1 2 4])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(3))
  [h(1),labs(3,:)]=plott2(data);
  [x,y,z]=latitude(90-th0,1);
  ha(2)=plot3(x,y,z,'y-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(4))
  [h(2),labs(4,:)]=plott2(data2);
  hold on
  [x,y,z]=sph2cart(-gam/180*pi,(90-bta)/180*pi,1.01);
  ps(1)=plot3(x,y,z,'s');
  [x,y,z]=latitude(bta+90,1);
  h(16)=plot3(x,y,z,'w-');
  [x,y,z]=longitude(-gam,1);
  h(11)=plot3(x,y,z,'w-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(5))
  [h(3),h(4)]=plott3(data);
  hold on
  [x,y,z]=latitude(90-th0,3);
  ha(3)=plot(x,y,'y-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(6))
  [h(5),h(6)]=plott3(data2);
  hold on
  [xgb,ygb]=mollweide(-gam*pi/180,(90-bta)*pi/180,pi);
  ps(4)=plot(xgb,ygb,'s');
  [x,y]=longitude(-gam,3);
  h(12)=plot(x,y,'w');
  [x,y]=latitude(bta+90,3);
  h(13)=plot(x,y,'w');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(7))
  [h(7),h(8)]=plott4(data);
  hold on
  [x,y,z]=latitude(90-th0,2);
  ha(1)=plot(x,y,'y-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(8))
  [h(9),h(10)]=plott4(data2,[alp bta gam]);
  hold on
  ps(5)=plot(-gam,bta+90,'s');
  [x,y,z]=latitude(bta+90,2);
  h(14)=plot(x,y,'w-');
  [x,y,z]=longitude(-gam,2);
  h(15)=plot(x,y,'w-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  longticks(ah(7:8))
  undeggies(ah(7:8))
  set(ah(7:8),'xtick',[0:90:360],'xtickl',[0:90:360],'xlim',[0 360],...
	      'ytick',[-90:45:90],'ytickl',[-90:45:90],'ylim',[-90 90])
  deggies(ah(7:8))
  set(ah(1:2),'camerav',8.5)
  set(ah(3:4),'camerav',8.5)
  set(ah(5:6),'camerav',5)
  set(h,'LineW',0.5)
  set(ps([1 4 5]),'MarkerF','g','MarkerE','k')
  t(1)=supertit(ah(1:2),sprintf('%s= %i ; %s= %i ; %s= %i',...
				'\alpha',alp,'\beta',bta,'\gamma',gam));
  t(2)=title(ah(8),sprintf('latitude %i ; longitude %i',...
			   bta+90,-gam));
  fig2print(gcf,'tall')
  figdisp([],demonr)
  set(gcf,'inv','off','color','w')
end

% Auxiliary plotting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spat,xyz]=plott1(cofs)
spat=plotplm(cofs,[],[],2,1); axis image; hold off
plotonsphere(spat,0.2); axis image on; 
xyz(1)=xlabel('x');
xyz(2)=ylabel('y');
xyz(3)=zlabel('z'); grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ch,xyz]=plott2(spat)
ch=plotonearth(spat,1); axis image on; box on ; grid on
set(ch,'Color','w')
xyz(1)=xlabel('x');
xyz(2)=ylabel('y');
xyz(3)=zlabel('z'); grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ch,ph]=plott3(spat)
[jk,chh,ph]=plotplm(spat,[],[],1,1); axis image; hold off
ch=chh{1}; set([ch ph],'Color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ch,ph]=plott4(spat,abg)
imagef([0 90],[360 -90],spat)
[jk,ch,lola]=plotcont;
[ph,lolap]=plotplates;
if nargin==2
  delete(ch); delete(ph)
  abg=abg*pi/180;
  % See if ROTTP behaves consistently with PLM2ROT
  [thetap,phip]=rottp(pi/2-lola(:,2)*pi/180,lola(:,1)*pi/180,...
		      abg(1),abg(2),abg(3));
  [thetapp,phipp]=rottp(pi/2-lolap(:,2)*pi/180,lolap(:,1)*pi/180,...
		      abg(1),abg(2),abg(3));
  % Get rid of the unsightly connecting lines
  phip=[phip+(phip<0)*2*pi]*180/pi;
  thetap=90-thetap*180/pi;
  [thetap,phip]=penlift(thetap,phip);
  phipp=[phipp+(phipp<0)*2*pi]*180/pi;
  thetapp=90-thetapp*180/pi;
  [thetapp,phipp]=penlift(thetapp,phipp);
  hold on
  ch=plot(phip,thetap,'b');
  ph=plot(phipp,thetapp,'b');
  hold off
end
set([ch ph],'Color','k')
