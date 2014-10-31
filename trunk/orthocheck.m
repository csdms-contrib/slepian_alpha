function [ngl1,ngl2,com,Vc,nofa,zmean]=orthocheck(C,V,TH,m,sord,ntw,cmean)
% [ngl1,ngl2,com,Vc,nofa,zmean]=ORTHOCHECK(C,V,TH,m,sord,ntw,cmean)
%
% Checks the orthonormality of a SINGLE-ORDER spherical harmonic
% expansion over the UNIT SPHERE and over a SINGLE or DOUBLE spherical CAP.
% If no EIGENVALUES are known, calculates them by integration.
% All in fully unit-normalized real harmonics (to 1 over unit sphere).
%
% INPUT:
%
% C        Spherical harmonic expansion coefficients (sine or cosine)
% V        Eigenvalues, in a matrix, if known
% TH       Radius of the spherical polar cap, in RADIANS
% m        Angular order, m>=0
% sord     1 Single polar cap [default]
%          2 Double polar cap
% ntw      0 Worry about orthogonality and report deviance [default]
%          1 Don't and simply calculate eigenvalues
% cmean    Also check the mean over the entire sphere [default: 0]
%
% OUPUT:
%
% ngl1     Number of GL points on the unit sphere to satisy orthogonality
% ngl2     Number of GL points on the domain to satisy orthogonality
% com      Center of mass of SINGLE cap (see Freeden/Narcowich)
% Vc       Eigenvalues deduced from the integration
% nofa     Normalization factor to make eigenfunctions integrate to unity
% zmean    The mean of the eigenfunction over the entire sphere
%
% Displays a message with the mean absolute error over the coefficients,
% the unit sphere, and the domain. The number of GL points are obtained
% by assuming that the orthogonality is satisfied and then determining the
% number which best verifies this, which is useful for non-zero order.
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

defval('sord',1)
defval('ntw',0)
defval('cmean',0)

% Set tolerance level
tol=10^-12;

% Figure out the maximum degree of the expansion
L=size(C,1)+abs(m)-1;

% Verify normalization of the coefficients; this should be zero
% This here is testing SDW (2006) equation (4.8), first part
% This here is testing SDW (2006) equation (5.8), first part
err(1)=mean(mean(abs(C'*C-eye(size(C,2)))));
% Hmm... why not also test C'*K*C'? Should be the eigenvalues.

% Verify normalization to unity over the unit sphere
% This here is testing SDW (2006) equation (4.8), first part
% Upping the number of GL points is crucial
if m==0
  % At least two integration points, to cover the L=0 case without messing
  % with plm2th
  ngl=max(2*L,2);
  % Normalization over the unit sphere
  [w,x]=gausslegendrecof(ngl,[],[-1 1]);
  % This is SDW (2006) equation (5.4) at the GL points
  % BUT IT ALSO INCLUDES sqrt(2-dom) of equation (5.9)
  EGL=pl2th(C,acos(x),1);
  % This second error should be zero
  [err(2),orv2,zmean]=orvo1(w,EGL,m,cmean);
  % If it is, don't worry about orv2, which is the identity, if not, keep it
  if err(2)<=tol ; orv2=NaN ; end
  
  ngl1=ngl;
  if sord==1
    % Center of mass
    com=2*pi*diag(EGL'*diag(w.*x(:))*EGL);
    com=com';
  else
    com=0;
  end
elseif m~=0
  if ntw==0
    ngl=2*L+[1:2:100];
  else
    ngl=2*L+1;
  end

% Iteratively determines the best number of integration points
  for index=1:length(ngl)
    [w,x]=gausslegendrecof(ngl(index),[],[-1 1]);
    % This is SDW (2006) equation (5.4) at the GL points
    % BUT IT ALSO INCLUDES sqrt(2-dom) of equation (5.9)
    EGL=plm2th(C,acos(x),m,1);
    [trerr(index),orv2,zmean]=orvo1(w,EGL,m,cmean);
    if trerr(index)<tol; break; end
  end
  [err(2),jk]=min(trerr);
  % If it is, don't worry about orv2, which is the identity, if not, keep it
  if err(2)<=tol ; orv2=NaN ; end

  ngl1=ngl(index);
  com=NaN;
  if index==length(ngl) & index>1; 
    disp('Could do with more integration accuracy')
  end
end

clear trerr

% Verify normalization to eigenvalue over the polar cap(s)
% This here is testing SDW (2006) equation (4.8), second part
% If ngl is too big leads to inaccurate results!
% If V is a matrix, it's the standard case, and we check.
% If V is a vector then it's the Sturm-Liouville eigenvalue and we
% need to still calculate the localization eigenvalue.
intv=[cos(TH) 1];
if m==0
  [w,x]=gausslegendrecof(ngl,[],intv);
  % This is SDW (2006) equation (5.7) at the GL points
  % BUT IT DOES INCLUDE sqrt(2-dom) of equation (5.9)
  EGL=pl2th(C,acos(x),1);
  if sord==1 % SIMPLE CAP
    [err(3),Vc]=orvo2(w,EGL,m,V,orv2);
  elseif sord==2 % DOUBLE CAP
    % The top half
    % The sqrt(2) had already been taken care of, thus the longitudinal
    % integral is either pi or 2pi
    orv4=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
    % Add the bottom half, same weight, flipped x
    EGL=pl2th(C,acos(-flipud(x)),1);
    % The sqrt(2) had already been taken care of, thus the longitudinal
    % integral is either pi or 2pi
    orv4=orv4+[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
    if prod(size(V))~=length(V) & ~isempty(V)
      % Then, not from GRUNBAUM2, eigenvalues known
      err(3)=mean(mean(abs(orv4-V)));
      Vc=V;
    else
      % Add it in for GRUNBAUM2, calculate them yourself
      if isnan(orv2) % ORTHOGONAL SYSTEM
	Vc=diag(orv4)';
	% But the off-diagonal terms must be zero for orthogonality
	err(3)=mean(mean(abs(orv4-diag(Vc))));
      else % NON-ORTHOGONAL SYSTEM
	Vc=diag(orv4./orv2)';
	err(3)=1;
      end
    end
  end
  ngl2=ngl;
elseif m~=0
  % Iteratively determines the best number of integration points
  for index=1:length(ngl)
    [w,x]=gausslegendrecof(ngl(index),[],intv);
    % This is SDW (2006) equation (5.7) at the GL points
    % BUT IT DOES INCLUDE sqrt(2-dom) of equation (5.9)  
    EGL=plm2th(C,acos(x),m,1);
    if sord==1 % SIMPLE CAP
      [trerr(index),Vc]=orvo2(w,EGL,m,V,orv2);
      if trerr(index)<tol; break; end
    elseif sord==2 % DOUBLE CAP
      % The top half
      % The sqrt(2) had already been taken care of, thus the longitudinal
      % integral is either pi or 2pi
      orv4=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
      % Add the bottom half, same weight, flipped x
      EGL=plm2th(C,acos(-flipud(x)),m,1);
      % The sqrt(2) had already been taken care of, thus the longitudinal
      % integral is either pi or 2pi
      orv4=orv4+[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
      if prod(size(V))~=length(V) & ~isempty(V)
	% Then, not from GRUNBAUM2, eigenvalues known
	trerr(index)=mean(mean(abs(orv4-V)));
	Vc=V;
	if trerr(index)<tol; break; end
      else
	% ORTHOGONAL SYSTEM
	if isnan(orv2)
	  % Add it in for GRUNBAUM2
	  Vc=diag(orv4)';
	  % But the off-diagonal terms must be zero for orthogonality
	  trerr(index)=mean(mean(abs(orv4-diag(Vc))));
	  if trerr(index)<tol; break; end
	else % NON-ORTHOGONAL SYSTEM
	  Vc=diag(orv4./orv2)';
	  trerr(index)=1;
	  % No way of increasing accuracy here since result is not known
	  break
	end
      end
    end
  end
  [err(3),jk]=min(trerr);
  ngl2=ngl(index);
  if index==length(ngl) & index>1; 
    disp('Could do with more integration accuracy')
  end
end
% Report on the errors if you're interested only
if ntw==0
  if any(err>tol)
    warning(sprintf(...
	'ORTHOCHECK Normalization criteria NOT satisfied; mean errror %8.3e',...
	mean(err)))
  else
    disp(sprintf(...
	'ORTHOCHECK Normalization criteria satisfied to %8.3e',...
	mean(err)))
  end
  nofa=NaN;
else
  nofa=sqrt(diag(orv2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [erro,orv2,zmean]=orvo1(w,EGL,m,cmean)
% Normalization (and mean) over the unit sphere
% Calculate the latitudinal and longitudinal integral
% For the diagonal elements only?
% orv1=w(:)'*[EGL.^2]*(1+(m==0))*pi;
% For the entire matrix
% The sqrt(2) had already been taken care of, thus the longitudinal
% integral is either pi or 2pi
orv2=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
erro=mean(mean(abs(orv2-eye(size(orv2)))));
% Additionally check the mean of the function which should be zero in the
% bandpass case. But this means the mean of the radial function should be
% zero only for m=0 as for any other m the {cos,sin}(m*phi) takes care
if cmean==1
  zmean=w'*EGL*(m==0);
else
  zmean=NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [erro,Vc]=orvo2(w,EGL,m,V,orv2)
% Normalization over the domain
% Calculate the latitudinal and longitudinal integral
% For the diagonal elements only?
% orv3=w(:)'*[EGL.^2]*(1+(m==0))*pi;
% For the entire matrix
% The sqrt(2) had already been taken care of, thus the longitudinal
% integral is either pi or 2pi
orv4=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
if prod(size(V))~=length(V) & ~isempty(V)
  erro=mean(mean(abs(orv4-V)));
  Vc=V;
else
  % We don't check but calculate the eigenvalues as the diagonal
  if isnan(orv2)
    Vc=diag(orv4)';
    % But the off-diagonal terms must be zero for orthogonality
    erro=mean(mean(abs(orv4-diag(Vc))));
  else
    Vc=diag(orv4./orv2)';
    erro=0;
  end
end


