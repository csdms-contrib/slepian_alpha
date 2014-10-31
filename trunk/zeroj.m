function [s,C,S,L]=zeroj(l1,l2,l3,L,meth,C,S)
% [s,C,S,L]=ZEROJ(l1,l2,l3,L,meth,C,S)
%
% Wigner 3j-symbol from a database precomputed by WIGNERCYCLE...
% with bottom row (spherical harmonic orders) all zeros!
%
% INPUT:
%
% l1,l2,l3     Top row of the Wigner 3j symbol [may be vectors]
% L            The bandwidth of the database [default: best available]
% meth         0 Uses sparse matrices [elegant, but slow]
%              1 Performs linear search on unsorted array [slow]
%              2 Performs binary search on presorted array [default]
% C,S          The column and element vectors resulting from a previous load
%
% OUTPUT:
%
% s            The (vector of) Wigner 3j symbols with zero bottom row
% C,S          The column and element vectors good for the next load
% L            The L of the database that is loaded, should you not know
%
% EXAMPLE:
%
% zeroj('demo1') % Should return nothing if it all works
% [jk,C0,S0]=zeroj(0,0,22) % Loads best available database for later use
% [jk,C0,S0]=zeroj(0,0,0,22) % Loads/Makes the database for later use
%
% l=32; p=56;
% difer(1-sum((2*[0:l+p]+1).*zeroj(l,p,[0:l+p]).^2))
%
% SEE ALSO: THREEJ, SIXJ, WIGNERCYCLE
%
% Last modified by fjsimons-at-alum.mit.edu, 04/09/2007

if ~isstr(l1) % Not a demo
  defval('C',[])
  defval('S',[])
  
  % Method
  defval('meth',2)
  % disp(sprintf('Using method %i',meth))
  
  % All saved values must be integers
  if sum(mod(l1,1)) || sum(mod(l2,1)) || sum(mod(l3,1))
    error('All degrees must be integers for the database')
  end
  
  if isempty(C) && isempty(S)
    % Got something to do
    % What is the lowest of the available database bandwidths?
    Els=ls2cell(fullfile(getenv('IFILES'),'WIGNER','WIGNER0JCS-*-C'));
    % Here it breaks if there is NOTHING
    for index=1:length(Els)
      EL(index)=str2num(rindeks(parse(Els{index},'-'),2));
    end
    EL=sort(EL);
    % Bandwidth of the database; keep this as low as possible
    % Need to provide for empties if not found
    fmax=find(max([l1(:)' l2(:)' l3(:)'])<=EL);
    if ~isempty(fmax)
      defval('L',EL(indeks(fmax,1)))
    else
      defval('L',-1)
      % Maybe here need to do wignercycle and have a go again?
    end
    % Check, identify and load database
    if any([l1(:)' l2(:)' l3(:)']>L)
      error('Insufficient bandwidth for database')
      % HERE WE SHOULD RUN WIGNERCYCLE
    end
    fnpl1=sprintf('%s/WIGNER0JCS-%i-C',...
		  fullfile(getenv('IFILES'),'WIGNER'),L);
    fnpl2=sprintf('%s/WIGNER0JCS-%i-S',...
		  fullfile(getenv('IFILES'),'WIGNER'),L);
    if exist(fnpl1,'file')==2 && exist(fnpl2,'file')==2
      fmt1='uint64'; fmt2='float64';
      disp(sprintf('Loading %s',fnpl1))
      C=eval(sprintf('loadb(''%s'',''%s'')',fnpl1,fmt1));
      disp(sprintf('Loading %s',fnpl2))
      S=eval(sprintf('loadb(''%s'',''%s'')',fnpl2,fmt2));
      
      % Whatever happens, this better be sorted; check some entries
      randin=unique(ceil(rand(min(100,length(C)),1)*length(C)));
      if ~all(unique(C(randin))==C(randin))
	disp('Column arrays were not properly sorted')
	[C,j]=sort(C);
	writeb(C,fnpl1,fmt1)
	S=S(j);
	writeb(S,fnpl2,fmt2); clear j
	disp('Column arrays now sorted once and for all')
      end
    else
      % Precompute the database - with bottom rows of zero
      wignercycle(L,0);
      % And have a go again
      s=zeroj(l1,l2,l3,L);
    end
  else
    if isempty(L)
      error('If supplying vectors with coefficients must also specify bandwidth')
    end
    % Else have C and S from a previous load and do nothing, but check
    if any([l1(:)' l2(:)' l3(:)']>L)
      error('Insufficient bandwidth for database')
    end
  end

  % Find running index into this matrix
  % So you can make a vector of l1 l2 l3 and have all zero m's
  rind=addmabout(L,l1)+(L+1)^2*(addmabout(L,l2)-1)+...
       (L+1)^4*(addmabout(L,l3)-1);

  % Initialize results vector
  s=repmat(NaN,1,length(rind));

  switch meth
    case 0
     % Turn it into the properly indexed sparse matrix
     % This step takes some of time but most of all, memory
     W=sparse(1,C,S,1,(L+1)^8);

     % Extract the Wigner 3j-symbol
     s=W(1,rind);
   case 1
    for index=1:length(rind)
      posi=find(C==rind(index));
      if ~isempty(posi)
	s(index)=S(posi);
      else 
	s(index)=0;
      end
    end
   case 2
    % Binary search algorithm on sorted arrays
    for index=1:length(rind)
      posi=binsearch(C,rind(index));
      if ~isempty(posi)
	s(index)=S(posi);
      else 
	s(index)=0;
      end
    end
   otherwise
    error('Specify valid method')
  end
elseif strcmp(l1,'demo1')
   difer(threej(16,16,6,0,0,0)-zeroj(16,16,6))
   difer(threej([2 0],[1 1],[1 1],0,0,0)-zeroj([2 0],[1 1],[1 1]))
   difer(threej(0,1,1,0,0,0)-zeroj(0,1,1))
   difer(threej(2,1,1,0,0,0)-zeroj(2,1,1))
   difer(threej(20,10,10,0,0,0)-zeroj(20,10,10))
   difer(threej(20,8,12,0,0,0)-zeroj(20,8,12))
   difer(threej(30,21,12,0,0,0)-zeroj(30,21,12))
   difer(threej(27,19,8,0,0,0)-zeroj(27,19,8))
   difer(threej(0:40,10,13,0,0,0)-zeroj(0:40,10,13))
end

