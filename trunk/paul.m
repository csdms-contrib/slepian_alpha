function varargout=paul(Lmax,x0)
% [Itab,Itabo,dels,dems]=paul(Lmax,x0)
%
% Recursively evaluates a list of integrals of a single Schmidt
% semi-normalized real Legendre polynomial P_lm(x)dx from x0 to 1.
% The iteration is set up according to Paul (1978) 'Recurrence relations
% for integrals of associated Legendre functions'.
%
% INPUT:
%
% Lmax          Bandwidth (maximum spherical harmonic degree)
% x0            Vector of lower bounds for the integrals 
%
% OUTPUT:
% 
% Itab          List of integrals of LEGENDRE(l,x,'sch') from x0 to 1 
% Itabo         Same but without the sqrt(2-dom) factor subsumed by Matlab
%
% EXAMPLE: 
%
% paul('demo1') Random Lmax, 100 values for x0, PAUL vs GAUSSLEGENDRE 
% paul('demo2') Complete list comparison PAUL vs GAUSSLEGENDRE
%
% Last modified by plattner-at-princeton.edu, 05/17/2011
% Last modified by fjsimons-at-alum.mit.edu, 06/01/2011

defval('Lmax',ceil(rand*256))
defval('x0',rand(256,1)*2-1)

if ~isstr(Lmax)

  % Turn the lower bounds into row
  x0=x0(:)';
  
  % first we test if Lmax is zero or positive
  if(Lmax<0)
    error('Lmax must be non-negative');
    Itable=[];
  end
      
  % Initialize Itab
  [Itab,Ptab]=deal(zeros(Lmax*(Lmax+1)/2+Lmax+1,length(x0)));
  
  % Now the starting evaluations for the LEGENDRE(l,x,'sch')
  Itab(1,:)=legendreint01(0,0,x0);
  Itab(2,:)=legendreint01(1,0,x0);
  % Get rid of the sqrt(2) factor
  Itab(3,:)=legendreint01(1,1,x0)/sqrt(2);

  % Now the Legendre functions themselves
  Ptab(2:3,:)=legendre(1,x0,'sch');
  % Get rid of the sqrt(2) factor
  Ptab(3,:)=Ptab(3,:)/sqrt(2);
  
    % This is a common argument
  y2=(1-x0.^2);
  
  % Now begin with the iterations 
  for l=2:Lmax
    % This is the index for m=0 at this l, i.e. addmup(l-1)+1
    lind=l*(l+1)/2+1;
    % This is the index for m=0 at l-1, i.e. addmup(l-2)+1
    lond=(l-1)*l/2+1;
    % This is the index for m=0 at l-2, i.e. addmup(l-3)+1
    lund=(l-2)*(l-1)/2+1;
    
    % Calculate the associated Legendre functions for this l
    Ptab(lind:lind+l,:)=legendre(l,x0,'sch');
    % Get rid of the sqrt(2) factor
    Ptab(lind+1:lind+l,:)=Ptab(lind+1:lind+l,:)/sqrt(2);
    
    % Some common factors at this l
    fux=(2*l-1)/(l+1);
    fex=sqrt((2*l-1)*(2*l-2)*2*l);

    % These are the cases where m~=l
    for m=0:l-1
      % Some common factors at this m
      fax=sqrt((l  +m)*(l  -m));
      fox=sqrt((l-1+m)*(l-1-m));
      if m<l-1
        Itab(lind+m,:)=Itab(lund+m,:)*(l-2)/(l+1)/fax*fox;
      end
      Itab(lind+m,:)=Itab(lind+m,:)+y2.*Ptab(lond+m,:)/fax*fux;
    end          
    
    % This is the case where m=l
    Itab(lind+l,:)=Itab(lund+l-2,:)*l*sqrt(2*l-3)/fex*fux-...
        y2.*Ptab(lond+l-2,:)/fex*fux;
  end
  
  [dems,dels,mz]=addmon(Lmax);

  % Now put the sqrt(2) factor in again since that's what we asked for
  nonzonal=skip(1:length(dems),mz);
  Itabo=Itab;
  Itab(nonzonal,:)=Itab(nonzonal,:)*sqrt(2);
  
  % Optional output
  varns={Itab,Itabo,dels,dems};
  varargout=varns(1:nargout);
elseif strcmp(Lmax,'demo1')
  % Demo to compare the Paul iteration method to Gauss-Legendre integration
  
  L=round(256*rand);
  m=round(L*rand);
  
  N=100;
  disp(sprintf('Calculate integrals for l up to %i for %i values for x0',L,N))
  x0=linspace(-1,1,N);
  
  % Paul's (1978) method 
  tic
  Itab=paul(L,x0);
  toc    
  
  % Gauss Legendre integration of the associated Legendre polynomials       
  tic
  %integrand=inline(sprintf(...
  %    ['rindeks(legendre(%i,x,''sch''),%i)'],...	
  %    L,m+1));
  integrand=@(x) rindeks(legendre(L,x,'sch'),m+1);
  %[w,xgl,nsel]=gausslegendrecof(max(L,200*((m)~=0)));
  [w,xgl,nsel]=gausslegendrecof(1000);
  in=zeros(1,N);
  for i=1:N
    in(i)=gausslegendre([x0(i) 1],integrand,[w(:) xgl(:)]);
  end
  toc
  
  subplot(211)
  plot(x0,Itab(L*(L+1)/2+m+1,:),'r-+',x0,in,'b')
  xlabel('Lower limit x0')
  ylabel('I_{l,m}')
  title(sprintf('l = %d, m = %d',L,m))
  legend('paul','gausslegendre')
  
  % Calculate the error
  errs=Itab(L*(L+1)/2+m+1,:)-in;
  difer(errs);
  
  subplot(212)
  plot(x0,errs,'k')
  xlabel('Lower limit x0')
  ylabel('Difference')
  title(sprintf('l = %d, m = %d',L,m))
elseif strcmp(Lmax,'demo2') 
  
  L=round(rand*50);
  
  N=10;
  x0=-1:2/(N-1):1;
  x0=x0(:)';
  
  % Paul's method 
  tic
  Itab=paul(L,x0);
  toc

  % Gauss Legendre integration of the associated Legendre polynomials     
  gleg=zeros(L*(L+1)/2+L+1,N);
  [w,xgl,nsel]=gausslegendrecof(2000);
  tic
  for LL=0:L
    for m=0:LL
      % Slower...
      %integrand=inline(sprintf(...
      %    ['rindeks(legendre(%i,x,''sch''),%i)'],...	
      %    LL,m+1));
      % Faster...
      integrand=@(x) rindeks(legendre(LL,x,'sch'),m+1);
      in=zeros(1,N);
      for i=1:N
        in(i)=gausslegendre([x0(i) 1],integrand,[w(:) xgl(:)]);
      end
      gleg(LL*(LL+1)/2+m+1,:)=in;
    end
  end
  toc
  
  % Calculate errors of everything
  err=Itab-gleg;
  difer(err,7)

  % Setup error matrix
  errs=zeros(L,L);
  for LL=0:L
    for m=0:LL
      errs(LL+1,m+1)=sum(abs(gleg(LL*(LL+1)/2+m+1,:)...
                             -Itab(LL*(LL+1)/2+m+1,:)));
    end
  end

  clf
  imagefnan([0 0],[L L],errs,'kelicol',minmax(errs))
  longticks(gca,2)
  xlabel('order m'); ylabel('degree l');
  [cb,xcb]=addcb('hor',minmax(errs),minmax(errs),'kelicol',max(errs(:))); 
  movev(cb,-.05)
  set(xcb,'string',sprintf('sum of the absolute errors for %i lower limits',N))
end
