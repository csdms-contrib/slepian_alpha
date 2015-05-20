function varargout=legendreint01(l,m,x0)
% val=legendreint01(l,m,x0)
%
% Evaluates integrals of a Schmidt semi-normalized real Legendre
% polynomial P_lm(x)dx between x0 to 1 for l=0 or 1, analytically.
%
% INPUT:
% 
% l         Angular degree of the polynomial, l>=0
% m         Angular order of the polynomial, 0<=m<=l
% x0        Lower bound(s) for the integral(s)
%
% OUTPUT: 
%
% val       Integral(s) of LEGENDRE(N,x,'sch') or LIBBRECHT(N,x,'sch')
%
% SEE ALSO:
%
% LIBBRECHT, PAUL, LEGENREPRODINT
%
% EXAMPLE:
%
% legendreint01('demo1') Compare analytic to numerical Legendre functions
% legendreint01('demo2') Compare Gauss-Legendre integration to analytic
% legendreint01('demo3') Evaluates Gauss-Legendre integration errors 
%
% Last modified by plattner-at-princeton.edu, 05/17/2011
% Last modified by fjsimons-at-alum.mit.edu, 06/01/2011

% Calculates the integral of a single Legendre function between x0 and 1
if ~isstr(l)
  
  if l<0 || m<0 || m>l
    error('Bad choice of l and m');
  end
  
  if m>1
    error('Use PAUL or LEGENDREPRODINT instead.')
  end

if l==0
  % l=0 m=0, the integral of 1
  val=(1-x0);
elseif l==1
  if m==0
    % l=1 m=0, the integral of x
    val=(1-x0.^2)/2; 
  else
    % l=1 m=1, the integral of sqrt(1-x.^2)
    % See Research Book 9 page 79 for these alternative expressions
    %val=(   +acos(x0)/2-(x0.*sqrt(1-x0.^2))/2);
    %val= coscos([],1,1,[asin(x0(:)) repmat(pi/2,length(x0),1)])';
    %val=-sinsin([],1,1,[acos(x0(:)) repmat(0,length(x0),1)])';
    val=(pi/4-asin(x0)/2-(x0.*sqrt(1-x0.^2))/2);
  end
end

% Optional output
varns={val};
varargout=varns(1:nargout);

elseif strcmp(l,'demo1')
    x=linspace(-1,1,100);

    clf
    subplot(311)
    plot(x,ones(size(x)),'k')
    hold on
    plot(x,legendre(0,x,'sch'),'o')
    plot(x,libbrecht(0,x,'sch'),'+')

    subplot(312)
    plot(x,x,'k')
    hold on
    plot(x,rindeks(legendre(1,x,'sch'),1),'o')
    plot(x,rindeks(libbrecht(1,x,'sch'),1),'+')

    subplot(313)
    plot(x,sqrt(1-x.^2),'k')
    hold on
    plot(x,rindeks(legendre(1,x,'sch'),2),'o')
    plot(x,rindeks(libbrecht(1,x,'sch'),2),'+')
elseif strcmp(l,'demo2')
    x=linspace(-1,1,100);
    ngl=5;
    
    integrand1=inline('ones(size(x))');
    integrand2=inline('x');
    integrand3=inline('sqrt(1-x.^2)');
    
    [gl1,gl2,gl3]=deal(zeros(size(x)));
    for i=1:length(x)
      % Here by GL integration
      gl1(i)=gausslegendre([x(i) 1],integrand1,ngl);
      gl2(i)=gausslegendre([x(i) 1],integrand2,ngl);
      gl3(i)=gausslegendre([x(i) 1],integrand3,ngl);
    end
    % Now the code itself
    val1=legendreint01(0,0,x);
    val2=legendreint01(1,0,x);
    val3=legendreint01(1,1,x);
    % And another piece of code at a random point
    randi=ceil(rand*length(x));
    opti={'automatic','dumb','gl','paul'};
    % Note that the choice may be overridden inside LEGENDREPRODINT
    rando=ceil(rand*length(opti));
    % Use various algorithms explicitly at random
    vol1=legendreprodint(0,0,0,0,x(randi),opti{rando});
    vol2=legendreprodint(1,0,0,0,x(randi),opti{rando});
    vol3=legendreprodint(1,1,0,0,x(randi),opti{rando});

    clf
    subplot(311)
    plot(x,val1,'k')
    hold on
    plot(x,gl1,'b+')
    plot(x(randi),vol1,'rv','MarkerF','r','markers',6)

    subplot(312)
    plot(x,val2,'k')
    hold on
    plot(x,gl2,'b+')
    plot(x(randi),vol2,'rv','MarkerF','r','markers',6)

    subplot(313)
    plot(x,val3,'k')
    hold on
    plot(x,gl3,'b+')
    plot(x(randi),vol3,'rv','MarkerF','r','markers',6)
elseif strcmp(l,'demo3')
    x=-1;
    ngl=linspace(1,1000,25);
    
    integrand1=inline('ones(size(x))');
    integrand2=inline('x');
    integrand3=inline('sqrt(1-x.^2)');
    
    [gl1,gl2,gl3]=deal(zeros(size(x)));
    for i=1:length(ngl)
      % Here by GL integration
      gl1(i)=gausslegendre([x 1],integrand1,ngl(i));
      gl2(i)=gausslegendre([x 1],integrand2,ngl(i));
      gl3(i)=gausslegendre([x 1],integrand3,ngl(i));
    end

    % Now the code itself
    val1=legendreint01(0,0,x);
    val2=legendreint01(1,0,x);
    val3=legendreint01(1,1,x);

    clf
    subplot(311)
    semilogy(ngl,abs(gl1-val1),'k')

    subplot(312)
    semilogy(ngl,abs(gl2-val2),'k')

    subplot(313)
    semilogy(ngl,abs(gl3-val3),'k')
end
