function W=wigner3j(l1,l2,l3,m1,m2,m3,method,recon)
% DO NOT USE THIS FUNCTION: SEE WIGNER0J, WIGNER3JM, THREEJ.
% MOST ANNOYINGLY, THE LAST ORDER SIGN IS UNINTUITIVE. ALSO, IT'S A
% MIXTURE OF WIGNER AND GAUNT COEFFICIENTS AND NORMALIZATIONS. SEE DEMO6.
% Used only for initialization of WIGNERJM/WIGNER0J.
%
% W=wigner3j(l1,l2,l3,m1,m2,m3,method,recon)
%
% Calculates absolute values of Wigner 3j coefficients in the context of
% UN-normalized associated Legendre functions. One half the inner product
% over the entire interval of three standard, unnormalized, associated
% Legendre polynomials of degree l1,l2,l3 and positive order m1,m2,m3
% is the square of this Wigner 3j symbol, Dahlen & Tromp Eq. (C.226).
% 
% INPUT:
%
% l1,l2,l3     Top row of the Wigner3j symbol; positive angular degree
% m1,m2,m3     Bottom row of the Wigner3j symbol
% method       'automatic' (default), or force:
%              'gl'    by Gauss-Legendre integration
%              'recu'  by recursion formula
%              Note that, depending on what's possible, you may not
%              always get what you asked for, it defaults back
% recon        1 Compute, maybe save to database [default]
%              0 Allow database reading
%
% OUTPUT:
%
% W            The result, perhaps with all the recursions as a vector.
%
% EXAMPLE:
%
% wigner3j('demo1') % Compare against LEGENDREPRODINT
% wigner3j('demo2') % Compare the two calculation methods
% wigner3j('demo3') % Test various normalizations
% wigner3j('demo4') % Test against an exact arithmetic calculator
% wigner3j('demo5') % Test the special formulas against WIGNER3JM
% wigner3j('demo6') % Response to a user query
%
% See also: WIGNER0J, WIGNER3JM
%
% Last modified by fjsimons-at-alum.mit.edu, 05/31/2011

defval('m1',0)
defval('m2',0)
defval('m3',0)
defval('method','automatic')
defval('recon',1)

% What's coming behind here more properly fits in WIGNER3JM and uses a
% different normalization than what is being used behind
if m1<0 || m2<0 || m3<0 
  % Only special cases are allowed, with the Condon-Shortley phase
  % Top row is of the form (l+a b l) with a=0,1,2, and b=1,2
  % Bottom row is of the form (-m 0 m)
  if l2==1 && m2==0 && m1==-m3 && l1==l3+1
    l=l3; m=abs(m1);
    W=(-1)^(l+m+1)*sqrt(2*((l+1)^2-m^2)/((2*l+3)*(2*l+2)*(2*l+1)));
    disp('WIGNER3J by special formula 1')
    return
  end
  if l2==2 && m2==0 && m1==-m3 && l1==l3+2
    l=l3; m=abs(m1);
    W=(-1)^(l+m)*sqrt(6*((l+1)^2-m^2)*((l+2)^2-m^2)/...
                        ((2*l+5)*(2*l+4)*(2*l+3)*(2*l+2)*(2*l+1)));
    disp('WIGNER3J by special formula 2')
    return
  end
  if l2==2 && m2==0 && m1==-m3 && l1==l3+1
    l=l3; m=abs(m1);
    W=(-1)^(l+m+1)*2*m*sqrt(6*((l+1)^2-m^2)/...
			    ((2*l+4)*(2*l+3)*(2*l+2)*(2*l+1)*2*l));
    disp('WIGNER3J by special formula 3')
    return
  end
  if l2==2 && m2==0 && m1==-m3 && l1==l3
    l=l3; m=abs(m1);
    W=(-1)^(l+m)*2*(3*m^2-l*(l+1))/...
                sqrt((2*l-1)*(2*l)*(2*l+1)*(2*l+2)*(2*l+3));
    disp('WIGNER3J by special formula 4')
    return
  else
    error('Better use WIGNER3JM, WIGNERCYCLE, and THREEJ ')
  end
end

% What's coming behind here has a very different normalization and no
% longer properly fits with the shortcuts inserted above. Now it treats
% the input orders as those of the triple Legendre product and constructs
% a Wigner-like symbol for it. Combined with the different normalization,
% this is not the best approach here, though it's still useful to check.
if ~isstr(l1)
  fnpl=sprintf('%s/W3J-%i-%i-%i-%i-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'WIGNER'),l1,l2,l3,m1,m2,m3);
  if exist(fnpl,'file')==2 && recon==0
    load(fnpl)
    %disp(sprintf('%s loaded by WIGNER3J',fnpl))
  else
    % Negative m's are allowed unless the result is computed by integration 
    if abs(m1)>l1 || abs(m2)>l2 || abs(m3)>l3
      error('Mismatch of angular degrees and orders')
    end 

    % Conditions of Dahlen and Tromp Eq. (C.186) % BUT NO NEGATIVE M
    if ~triangle(l1,l2,l3) || m1~=(m2+m3)
      W=0;
      return
    end

    % Conditions of Dahlen and Tromp Eq. (C.219)
    if m1==0 && m2==0 && m3==0 && mod(l1+l2+l3,2)
      W=0;
      return
    end    
    
    % Dahlen and Tromp Eqs. (B.49) and (C.226)
    % Effectively it is now a double, not a triple product
    if l3==0 && m1==m2 && ~strcmp(method,'recu') && ~strcmp(method,'gl') 
      % The triangle equality with l3=0 implies l1==l2
      % The initial check on m implies that m3==0
      if (abs(l1)+abs(m1))<21
	disp('WIGNER3J by normalization formula Dahlen and Tromp B.49')
	W=sqrt(factorial(l1+m1)/factorial(l1-m1)/(2*l1+1));
	return
      else
	% If factorial fails need to use Gauss-Legendre
	W=wigner3j(l1,l2,l3,m1,m2,m3,'gl');
	return
      end
    end

    % Via Gauss-Legendre integration, forced or by choice
    if (l1~=l2 || m1~=0 || m2~=0 || m3~=0 || strcmp(method,'gl')) ...
	    && ~strcmp(method,'recu')
      % If the m is odd then the associated is not a polynomial
      % and this is not the best way; would need to add new recursion
      % relations to cover this case.
      integrand=inline(sprintf(...
	  ['rindeks(legendre(%i,x),%i).*',...
	   'rindeks(legendre(%i,x),%i).*',...
	   'rindeks(legendre(%i,x),%i)'],...
	  l1,m1+1,l2,m2+1,l3,m3+1));
      % Sign information is lost, however.
      % I suppose you can put this back in from Luscombe
      % sign of last term is (-1)^(l2-l3+m2+m3), see WIGNER0J
      % No wait, it is DT C220
      % Dahlen and Tromp Eq. (C.226)
      W=sqrt(gausslegendre([-1 1],integrand,l1+l2+l3)/2);
      disp('WIGNER3J by Gauss-Legendre integration of triple product')
      return
    end

    % All-zero bottom row; Dahlen and Tromp Eq. (C.220); sign positive
    % Not by recursion
    if m1==0 && m2==0 && m3==0 && (l1+l2+l3)<20 && ~strcmp(method,'recu')
      S=(l1+l2+l3); 
      disp('WIGNER3J by summation formula Dahlen and Tromp C.220')
      % PUT THE SIGN IN THERE DUDE
      W=sqrt(1/factorial(S+1)*...
	     factorial(S-2*l1)*factorial(S-2*l2)*factorial(S-2*l3))*...
	factorial(S/2)/...
	factorial(S/2-l1)/factorial(S/2-l2)/factorial(S/2-l3);
      return
    end

    % All-zero bottom row and l1==l2; 
    % Recursion http://functions.wolfram.com/07.39.17.0013.01
    if m1==0 && m2==0 && m3==0 && l1==l2
      % By upward recursion and saving all the terms
      W=repmat(NaN,1,l3+1);
      % Dahlen and Tromp Eq. (B.49)
      W(1)=1/sqrt(2*l1+1); % Zeroth term
      for j=2:2:l3 
	% All odd terms zero from triangle rule with l1==l2 
	W(j)=0;
	% Keep sign positive
	W(j+1)=(j-1)/j*...
	       sqrt((4*l1^2+4*l1-j^2+2*j)/(4*l1^2+4*l1-j^2+1))...
	       *W(j-1);
      end
      disp('WIGNER3J by recursion formula Mathematica')
    end
  end
  
  % If you don't have anything by now, start all over with the defaults
  defval('W',wigner3j(l1,l2,l3,m1,m2,m3))
elseif strcmp(l1,'demo1')
  j=0; m=0;
  for l=0:100
    d(l+1)=wigner3j(l,l,j,m,m,m)^2-legendreprodint(l,m,l,m,-1,'gl')/2;
  end 
  clf
  plot(0:100,abs(d)/ep,'LineW',2)  
  xlabel('l')
  ylabel(sprintf('W_{ll}^0 difference %s','\times'))
  title('Difference between WIGNER3J and direct GL integration');
  grid on
elseif strcmp(l1,'demo2')
  m=0; l=round(rand*100);
  for j=0:100
    D=wigner3j(l,l,j,m,m,m);
    d(j+1)=D(end)-wigner3j(l,l,j,m,m,m,'gl');
    clear D
  end
  clf
  plot(0:100,abs(d)/eps,'LineW',2)  
  xlabel('j')
  ylabel(sprintf('W_{%i%i}^j difference %s eps',l,l,'\times'))
  title('Difference between WIGNER3J using recursion and GL integration');
  grid on
elseif strcmp(l1,'demo3')
  l=round(100*rand(1));
  disp(sprintf('Test normalization at degree %i',l))
  difer(wigner3j(l,l,0,0,0,0)^2*(2*l+1)-1,[],[],NaN)
  l=round(100*rand(1)); j=round(100*rand(1));
  disp(sprintf('Test normalization at degrees %i and %i',l,j))
  a=0; for i=0:l+j; a=a+wigner3j(l,j,i)^2*(2*i+1); end
  difer(a-1,[],[],NaN)
elseif strcmp(l1,'demo4')
  % Compare with the algorithm by Stone and Wood, which differs 
  % at most by the sign, e.g.
  disp('See http://www-stone.ch.cam.ac.uk/cgi-bin/wigner.cgi')
  difer(wigner3j(1,9,10)-sqrt(10/399),[],[],NaN)
  difer(wigner3j(5,0,5,0,0,0)-sqrt(1/11),[],[],NaN)
elseif strcmp(l1,'demo5')
  l=round(100*rand(1));
  disp(sprintf('Compare with WIGNER3JM at degree %i',l))
  difer(wigner3j(l+1,1,l,-1,0,1)-indeks(...
      wigner3jm(l+1,1,l,-1,0,1),'end'),[],[],NaN)
  l=round(100*rand(1));
  disp(sprintf('Compare with WIGNER3JM at degree %i',l))
  difer(wigner3j(l+2,2,l,-1,0,1)-indeks(...
      wigner3jm(l+2,2,l,-1,0,1),'end'),[],[],NaN)
  l=round(100*rand(1));
  disp(sprintf('Compare with WIGNER3JM at degree %i',l))
  difer(wigner3j(l+1,2,l,-1,0,1)-indeks(...
      wigner3jm(l+1,2,l,-1,0,1),'end'),[],[],NaN)
  l=round(100*rand(1));
  disp(sprintf('Compare with WIGNER3JM at degree %i',l))
  difer(wigner3j(l+0,2,l,-1,0,1)-indeks(...
      wigner3jm(l+0,2,l,-1,0,1),'end'),[],[],NaN)
elseif strcmp(l1,'demo6')
  % If you want something that is proportional to the triple product of 
  % spherical harmonics Y_{53},Y_{00},Y_{53}, you'll want to try
  %l=5; m=3;
  % But now let's do random values
  maxl=48;
  l=round(maxl*rand(1)); m=round(l*rand(1));
    disp(sprintf(...
	'Comparing with WIGNER3JM at degree %i and order %i %i',l,m))
  W1=wigner3j(l,0,l,m,0,m);
  % or else this which is just by a different method
  W2=wigner3j(l,l,0,m,m,0);
  % but then you also expect this to be proportional to 
  W3=indeks(wigner3jm(l,0,l,m,0,-m),'end');
  % and you should be able to also use
  W4=threej(l,0,l,m,0,-m);
  % and the relation between the first and the last two, which pairs
  % are identical, is due to an unfortunate mixture of normalization
  % conventions and a deliberate confusion of Wigner and Gaunt
  % coefficients. Bottom line is that this is why WIGNER3J is obsolete.
  % Note that also the sign in WIGNER3J is not properly calibrated
  format long
  abs([W1/sqrt(factorial(l+m)/factorial(l-m)) ;
   W2/sqrt(factorial(l+m)/factorial(l-m)) ; W3 ; W4])
  % But the conclusion is that the various procedures are completely
  % attuned to one another.
  % Lastly, check your integrals over the sphere, this should be one 
  W1=indeks(wigner0j(l,0,l),'end')...
      *indeks(wigner3jm(l,0,l,m,0,-m),'end');
  W2=indeks(wigner3jm(l,0,l,0,0,0),'end')...
      *indeks(wigner3jm(l,0,l,m,0,-m),'end');
  W3=zeroj(l,0,l)*threej(l,0,l,m,0,-m);
  (-1)^m*[W1 ; W2 ; W3]*sqrt((2*l+1)*(2*0+1)*(2*l+1))
  format short
end

