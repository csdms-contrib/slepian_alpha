function varargout=libbrecht(N,X,norma,tst,Mo)
% [P,dP]=libbrecht(N,x,norma,tst,Mo)
%
% Computes associated Legendre functions as a function of X=cos(theta)
% and their derivatives with respect to theta=acos(X). Compare LEGENDREDIFF.
% Note that the (-1)^m is not part of the definition of the option 'sch'
% in LEGENDRE, but that the 'sch' in Matlab contain a sqrt(2-dom). Since
% the derivative is with respect to acos(X) to get the derivative with
% respect to X we need to multiply this by -1/sin(acos(X)).
%
% INPUT:
%
% N           Scalar harmonic degree N>=0
% x           Vector argument in interval [-1 1]
% norma       'sch' SCHMIDT SEMI-NORMALIZED REAL:
%                   Pl0(1)=1 and norm for m==0 is 2/(2l+1).
%                   Plm(1)=0 and norm for m>0 is 4/(2l+1).
%                   cos/sin harmonics normalized to 4*pi/(2l+1)
%             'fnr' FULLY NORMALIZED REAL:
%                   Pl0(1)=sqrt(2*l+1) and norm for m==0 is 2 
%                   Plm(1)=0 and norm for m>=0 is 4 
%                   cos/sin harmonics normalized to 4*pi
%                   These 'sch'*sqrt(2l+1) are equivalent to 
%                   sqrt(2-dom)*sqrt(4*pi)*X*(-1)^m of XLM
%             'fnc' FULLY NORMALIZED COMPLEX [default]:
%                   P(1)=sqrt((2*l+1)/4/pi) and norm is 1/2/pi.
%                   Real spherical harmonics normalized to
%                   1 for m==0 and 1/2 for m>0
%                   Complex spherical harmonics normalized to 1.
%                   This then is identical to XLM.
%             'unn' UNNORMALIZED (not supported):
%                   Pl0(1)=1 and norm is 2/(2l+1)*(l+m)!/(l-m)!.
%                   cos/sin m==0 to 4*pi/(2l+1)*(l+m)!/(l-m)!
%                   cos/sin m>0  to 2*pi/(2l+1)*(l+m)!/(l-m)!
% tst         Perform simple normalization test using Simpson's rule
%             See also the exact integration of GAUSSLEGENDRE.
% Mo          A particular order M (but can't avoid computing all)
%
% OUTPUT
%
% P           The Legendre functions at the points X=cos(theta)
% dP          And their derivatives with respect to acos(X)=theta.
%
% Computes associated Legendre functions and their derivatives, of
% harmonic degree N and for all orders m=0, 1, ..., N, evaluated for each 
% element of X (where reals -1 <= X <= 1). N must be a scalar. The
% algorithm is reported to be stable up to N=500 or so.
%
% See Masters and Richards-Dinger, 1998 and Libbrecht 1985.
%
% EXAMPLE:
%
% [P,dP]=libbrecht(round(rand*20),linspace(-1,1,1000),'sch',[],1); plot(dP)
% libbrecht('demo1')    % Compares this algorithm with Schmidt LEGENDRE
% libbrecht('demo2')    % Plots Legendre function and its derivative
% libbrecht('demo3')    % Compares this algorithm with my LEGENDREDIFF
% libbrecht('demo4')    % Compares with asymptotic expression of BACKUS
% libbrecht('demo5')    % Compares with simple first difference
% libbrecht('demo6')    % Relation with XLM
% hilbxlm('demo1')      % Compares with asymptotic expression of HILB
% dahlen('demo1')       % Compares with asymptotic expression of DAHLEN
% legendrediff('demo1') % Compares with LEGENDREDIFF and first difference
%
% See also DAHLEN, HILBXLM, BACKUS, LEGENDREDIFF
%
% Last modified by fjsimons-at-alum.mit.edu, 09/08/2011

if ~isstr(N)
  defval('norma','fnc')
  defval('Mo',NaN)

  % Initialize arrays
  X=X(:)';
  if any(X>1 | X<-1)
    warning('X must contain real values between -1 <= X <= 1')
  end
  P=repmat(NaN,N+1,length(X));
  dP=repmat(NaN,N+1,length(X));
  sint=sin(acos(X));
  
  % Part of the normalization constant (see also below)
  switch norma
   case 'fnc'
    Kst1=(-1)^N*sqrt((2*N+1)/4/pi);
   case {'fnr'}
    Kst1=sqrt(2*N+1);
   case {'sch'}
    Kst1=1;
   case 'unn'
    error('Not supported (why should it)')
   otherwise
    error('Specify valid normalization scheme')
  end

  % Start pathological cases
  % Handle N=0
  if N==0
    P=repmat(Kst1,1,length(X));
    dP=repmat(0,1,length(X));
    % Optional output
    varns={P,dP};
    varargout=varns(1:nargout);
    return
  end

  % Here is the Libbrecht algorithm. Compute starting prefactor
  % sqrt(1/factorial(2*N))*factorial(2*N)/(2^N)/factorial(N)
  % by computing (1/2)*(3/4)*(5/6)*...*((2l-1)/(2l)):
  f1=1;
  for i=1:N
    f1=f1*(2*i-1)/(2*i);
  end
  f1=sqrt(f1);

  % Initial value for m=N, see MRD (1998) Eq. (6)
  P(N+1,:)=Kst1*f1;
  dP(N+1,:)=0;

  % Compute prefactors 
  M=1:N;
  f2=sqrt((N+M).*(N-M+1));

  % For all M downgoing (MRD (1998) Eq. (5))
  % Dividing out f2 progressively yields sqrt(1/factorial(2*N))
  % Note that you're switching sign here, too; which you need for the 
  % algorithm but need to undo for the Schmidt harmonics
  for m=N:-1:1
    P(m,:)=-(sint.*dP(m+1,:)+2*m*X.*P(m+1,:))/f2(m);
    dP(m,:)=sint.*P(m+1,:)*f2(m);
  end

  % Now convert back to ordinary spherical harmonics
  f3=1;
  for m=2:(N+1)
    dP(m,:)=(sint.*dP(m,:)+(m-1)*X.*P(m,:)).*f3;
    P(m,:)=P(m,:).*sint.*f3;
    f3=f3.*sint;    
  end

  % Conversion factors Pfnc*fac=Pxxx;
  switch norma
   case {'sch','fnr'}
    % Note that Matlab contains the sqrt(2-dom) factor as part of 'sch'
    xfax=sqrt(2);
    P(2:end,:)=P(2:end,:)*xfax;
    dP(2:end,:)=dP(2:end,:)*xfax;
    % Note that Matlab uses a definition for the associated Legendre
    % polynomials P(n,m;x) that includes the Condon-Shortley phase
    % convention (-1)^M, and that the Schmidt-normalized functions have
    % this factor yet again, so as to effectively get rid of it...
    % Undo the alternating sign here; last one was always positive
    for m=N:-2:1
      P(m,:)=-P(m,:);
      dP(m,:)=-dP(m,:);
    end
   case 'fnc'
    xfax=-1;
   otherwise
    error('Specify valid normalization')
  end

  % Handle very small arguments (at both poles)
  % Function oscillates wildly at the endpoints for high l
  % Fix this at the very end
  endpts=find(abs(sint)<eps); 
  if length(endpts)>2;
    error('Should not have more than two poles'); 
  end
  % Zeros where expected, see DT (B.64)
  dP(:,endpts)=0;
  P(2:end,endpts)=0;
  % Except at m==0  but watch to generate the correct sign;
  % P(1,X==1), the North pole, always needs to be positive by our
  % convention. If l is EVEN also the South pole is positive.
  % Somehow this logical construct popped out.
  for indx=1:length(endpts)
    % The condition which fixes the sign
    cnd=(-1)^(((2*sign(X(endpts(indx)))+2*(-1)^(N)+1)>=0)+1);
    % The Legendre functions at degree 0
    P(1,endpts(indx))=abs(Kst1)*cnd;
    % The derivatives of the Legendre functions at degree 1
    % dP(2,endpts(indx))=abs(Kst1)*sqrt(N*(N+1))/2*cnd;
    % Can't remember how I did this... check Bosch (2000)
    dP(2,endpts(indx))=abs(Kst1)*sqrt(N*(N+1))/2*xfax*cnd;
  end
  
  if ~isnan(Mo)
    P=P(Mo+1,:);
    dP=dP(Mo+1,:);
  end

  % Normalization test if it spans the entire interval
  defval('tst',[])
  if tst
    % This doesn't seem to work any longer - needs work!
    if prod(size(X))>1 && sum(abs(sort([X(1) X(end)])-[-1 1]))<eps
      f=(P.*P)'.*repmat(sint(:),1,N+1);
      ntest=simpson(acos(X),f);
      switch norma
       case 'fnr'
	shub=1/2/pi;
       case 'sch'
	shub=(4/(2*N+1));
	% Note that the SEMI-normalization refers to this;
	% The real cos/sin harmonics are all normalized to 4pi/(2l+1).
	ntest(1)=ntest(1)*2;
      end
      % Should integrate to shub but sign could be wrong depending on
      % direction of integration. Note that the m=0 term will never integrate
      % well using this dumb integration scheme.
      disp(sprintf('Average normalization error %8.6f %s',...
		   100*(1-sum(abs(ntest(1:end)/shub))/(N+1)),'%'))
    else
      disp('Could not verify normalization. See GAUSSLEGENDRE.')
    end
  end
  if any(isnan(P(:)) | isnan(dP(:)))
    error('Values returned are NaN - use Matlab''s LEGENDRE routine')
  end
  % Optional output
  varns={P,dP};
  varargout=varns(1:nargout);
elseif strcmp(N,'demo1')
  clf
  theta=linspace(-pi,pi,500);
  X=cos(theta); more off
  Ns=[0:255];
  err=nan(1,length(Ns));
  index=0;
  for N=Ns
    index=index+1;
    hsch=libbrecht(N,X,'sch');
    ksch=legendre(N,X,'sch');
    err(index)=sum(abs(hsch(:)-ksch(:)));
    %    plot(X,ksch,'k',X,hsch,'y'); title(num2str(N)); pause
  end
  plot(Ns,err,'LineW',2)
  shrink(gca,1.2,1.2)
  st=title('Associated Legendre functions (Schmidt): LIBBRECHT-LEGENDRE',...
	   'FontS',12);
  movev(st,5e-9)
  yl=ylabel('Maximum Absolute Error (m=0\rightarrowl)');
  xl=xlabel('Angular degree l');
  axis tight ; grid on
  longticks(gca,2)
  fig2print(gcf,'portrait'); id
elseif strcmp(N,'demo2')
  deg=round(rand*10);
  m=round(rand*deg);
  
  theta=linspace(0,pi,500);
  X=cos(theta); more off
  % If you give it real names it spits out output
  % Give different names to avoid varargout
  [P1,dP1]=libbrecht(deg,X,'sch');
  clf
  subplot(211)
  plot(acos(X),P1(m+1,:),'LineW',2)
  yl=ylabel('X_l^m(\theta)=N_l^m\timesP_l^m(cos\theta)');
  openup(gca,6)
  hold on
  yli=ylim;
  plot([pi/2 pi/2],yli,'k')
  sg={'odd','even'};
  title(sprintf(...
      'Associated Legendre functions (Schmidt); l=%i, m= %i, l+m %s',...
      deg,m,sg{2-mod(deg+m,2)}),'FontS',12)
  pilabels(gca)
  nolabels(gca,1); grid on; axis tight
  movev(gca,-.1)
  longticks(gca,2)
  xlim([-0.05 pi+0.05])
  subplot(212)
  plot(acos(X),dP1(m+1,:),'r','LineW',2)
  xl=xlabel('Colatitude (\theta)');
  yl=ylabel('dX_l^m(\theta)/d\theta');
  grid on; axis tight
  openup(gca,6)
  xlim([-0.05 pi+0.05])
  pilabels(gca)
  hold on
  yli=ylim;
  plot([pi/2 pi/2],yli,'k')
  longticks(gca,2)
  fig2print(gcf,'portrait'); id
elseif strcmp(N,'demo3')
  deg=round(rand*75);
  m=0;
  theta=linspace(0,pi,500);
  X=cos(theta); more off
  % If you give it real names it spits out output
  % Give different names to avoid varargout
  [P1,dP1]=libbrecht(deg,X,'sch');
  clf
  [PpL,PLm1,PL,PLp1]=legendrediff(deg,X);
  subplot(211)
  plot(acos(X),PL,'b','LineW',2)
  hold on
  plot(acos(X),P1(m+1,:),'y')
  yl=ylabel('X_l^m(\theta)=N_l^m\timesP_l^m(cos\theta)');
  movev(gca,-.1)
  legend('LEGENDREDIFF','LIBBRECHT')
  title(sprintf(...
      'Legendre functions (Schmidt); l=%i (m= %i)',...
      deg,m),'FontS',12)
  pilabels(gca)
  nolabels(gca,1); grid on; axis tight
  axis tight; openup(gca,6)
  longticks(gca,2)
  subplot(212)
  plot(acos(X),-sin(acos(X)).*PpL,'b-','LineW',2)
  hold on
  plot(acos(X),dP1(m+1,:),'y')
  xl=xlabel('Colatitude (\theta)');
  yl=ylabel('dX_l^m(\theta)/d\theta');
  pilabels(gca)
  axis tight; openup(gca,6)
  longticks(gca,2)
  grid on
  fig2print(gcf,'portrait'); id
elseif strcmp(N,'demo4')
  deg=round(rand*75); 
  m=0;
  theta=linspace(0,pi,500);
  X=cos(theta); more off
  % If you give it real names it spits out output
  % Give different names to avoid varargout
  P1=rindeks(libbrecht(deg,X,'sch'),1);
  yli=minmax(P1);
  [ap,th,th0,apb]=backus(deg,length(X),'sch');
  clf
  subplot(211)
  plot(acos(X),P1,'b','LineW',2)
  hold on
  plot(th,ap,'y')
  plot(th,apb,'r')
  pp(1)=plot([1 1]*th0,[-10 10],'k-');
  pp(2)=plot([pi pi]-th0,[-10 10],'k-');
  yl=ylabel('X_l^m(\theta)=N_l^m\timesP_l^m(cos\theta)');
  movev(gca,-.075)
  legend('LIBBRECHT','BACKUS','ROBIN')
  title(sprintf(...
      'Legendre functions (Schmidt); l=%i (m= %i)',...
      deg,m),'FontS',12)
  grid on; axis tight
  ylim(yli)
  pilabels(gca)
  longticks(gca,2)
  nolabels(gca,1); 
  subplot(212)
  plot(acos(X),abs(P1-ap),'y','LineW',2)
  hold on 
  plot(acos(X),abs(P1-apb),'r','LineW',2)
  pp(3)=plot([1 1]*th0,[-10 10],'k-');
  pp(4)=plot([pi pi]-th0,[-10 10],'k-');
  yl=ylabel('X_l^m(\theta)=N_l^m\timesP_l^m(cos\theta)');
  xl=xlabel('Colatitude (\theta)');
  axis tight
  pilabels(gca)
  ylim([-0.1 1]*1e-2) ; grid on
  set(pp,'Color',grey)
  longticks(gca,2)
  fig2print(gcf,'portrait'); id
elseif strcmp(N,'demo5')
  deg=round(rand*75);
  m=round(rand*deg);
  theta=linspace(0,pi,500);
  X=cos(theta); more off
  % If you give it real names it spits out output
  % Give different names to avoid varargout
  [P,dP]=libbrecht(deg,X,'sch',[],m);
  clf
  % The first one is zero anyway, so can ignore the warning
  warning off MATLAB:divideByZero
  plot(X,-dP./sin(theta),'b');
  warning on MATLAB:divideByZero
  hold on
  % Watch out since x is not equally spaced.
  plot(X(2:end)-indeks(diff(X),1)/2,...
       diff(P)./diff(X),'r')
  %ylim([-1 1])
  title(sprintf(...
      'derivative of Legendre function (Schmidt); l=%i (m= %i)',...
      deg,m),'FontS',12)
elseif strcmp(N,'demo6')
  N=100; l=round(rand*100);
  disp(sprintf('degree l = %i',l))
  P=libbrecht(l,linspace(-1,1,N),'fnr');
  X=xlm(l,[],acos(linspace(-1,1,N)))*sqrt(4*pi).*...
    repmat((-1).^[0:l]'.*sqrt(2-[0:l==0])',1,N);
  difer(P-X)
  P=libbrecht(l,linspace(-1,1,N),'fnc');
  X=xlm(l,[],acos(linspace(-1,1,N)));
  difer(P-X)
elseif strcmp(N,'demo7')
  % Products of derivatives of unequal degree but equal order when that
  % order is sign-switched
  mu=linspace(-1,1,100); 
  x7m3=xlm(7,-3,acos(mu));
  x8m3=xlm(8,-3,acos(mu));
  x73=xlm(7,3,acos(mu));
  x83=xlm(8,3,acos(mu));
  plot(mu(1:end-1),diff(x83).*diff(x73),'b')
  hold on
  plot(mu(1:end-1),diff(x8m3).*diff(x7m3),'r+')
  hold on
end
