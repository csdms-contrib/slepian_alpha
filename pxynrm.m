function norms=pxynrm(l,m,tol,fcn)
% norms=PXYNRM(l,m,tol,fcn)
%
% Checks the normalization of PLM, XLM, and YLM.
% DT (B.49/B.59)
%
% Not a stand-alone program. See PLM, XLM, and YLM.
%
% See Dahlen and Tromp (1998), Theoretical Global Seismology,
% DT (X.NN) refers to their numbered equations.

% Last modified by fjsimons-at-alum.mit.edu, 03/16/2016

norms=repmat(NaN,max(length(m),length(l)),max(length(m),length(l)));

% Normalization check at fixed order
if prod(size(m))==1
  % Check fixed-order normalization DT (B.49/B.59)
  for index=1:length(l)
    for ondex=index:length(l)
      % Check fixed-order normalization DT (B.49/B.59)
      [w,x]=gausslegendrecof(l(index)+l(ondex));
      switch fcn
       case 'P'
	PX1=plm(l(index),m,x);
	PX2=plm(l(ondex),m,x);
	fax=factorial(l+m)./factorial(l-m);
	adj=(ones+diag(1./(2./(2*l+1).*fax))-eye(size(norms)));
       case 'X'
	PX1=xlm(l(index),m,acos(x));
	PX2=xlm(l(ondex),m,acos(x));
	adj=(ones+diag(repmat(2*pi,1,length(l)))-eye(size(norms))); 
      end
      norms(index,ondex)=[PX1*diag(w)*PX2'];
      norms(ondex,index)=norms(index,ondex);
    end
  end
  norms=norms.*adj;
elseif prod(size(l))==1
  % Check fixed-degree same-order normalization
  [w,x]=gausslegendrecof(2*l);
  % Check fixed-degree different-order normalization, alas, with an
  % inaccurate quadrature rule
  for index=1:length(m)
    for ondex=index:length(m)
      mnz=[m(index) m(ondex)]; mnz=mnz(~~mnz);
      switch fcn
       case 'P'
	PX1=plm(l,m(index),x);
	PX2=plm(l,m(ondex),x);
	fax=factorial(l+m(index))./factorial(l-m(index));
	adj1=1./(2./(2*l+1).*fax);
	adj2=1./(2./(2*l+1));
	adj3=-(-1)^(m(index));
	adj4=(ones+diag(1/fax./mnz));
       case 'X'
	PX1=xlm(l,m(index),acos(x));
	PX2=xlm(l,m(ondex),acos(x));
	adj1=2*pi; 
	adj2=2*pi;
	adj3=0;
	adj4=(ones+diag(1./mnz))*(2*l+1)/4/pi;
      end	
      if abs(m(index))==abs(m(ondex))
	norms(index,ondex)=[PX1*diag(w)*PX2'];
	if m(index)==m(ondex)
	  % Diagonal term which should return 1
	  norms(index,ondex)=norms(index,ondex).*adj1;
	elseif m(index)==-m(ondex)
	  % Off-diagonal term but m=-m' thus returns 1 which we map to 0
	  norms(index,ondex)=norms(index,ondex).*adj2+adj3;
	end
      else % Off-diagonal term which should always return zero
	norms(index,ondex)=[PX1*diag(w./(1-x.^2))*PX2'];
	% Which is the nonzero one?
	norms(index,ondex).*adj4;
      end
      norms(ondex,index)=norms(index,ondex);
    end
  end
end
  
% Add it all up
chex=sum(sum(abs(norms-eye(size(norms)))));
if chex < tol
  disp(sprintf('PXYNRM: Normalization %0.2e < %0.2e satisfied',chex,tol))
else
  disp(sprintf('PXYNRM: Normalization %0.2e > %0.2e not satisfied',chex,tol))
end

