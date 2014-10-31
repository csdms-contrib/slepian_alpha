function varargout=plotoncube2(vname,dshell,tipe,kind,vscale,colmap,phil,setnan)
% [h,pg,pc]=PLOTONCUBE2(vname,dshell,tipe,kind,vscale,colmap,phil,setnan)
%
% Plots a SELECTED DEPTH of gridded whole-cubed-sphere data.
% A wrapper for PLOTONCUBE.
%
% INPUT:
%
% vname      A string identifying the vector [default: 'vecx.petsc'], or
%            A structure array with all of the information in there already
% dshell     An integer identifying the layer number, i.e. the depth shell
% tipe       See the option list in PLOTONCUBE
% kind       See the option list in PLOTONCUBE
% vscale     See the option list in PLOTONCUBE
% colmap     See the option list in PLOTONCUBE
% phil       See the option list in PLOTONCUBE
% setnan     See the option list in PLOTONCUBE
%
% OUTPUT:
%
% h          A series of six axis handles to the plots for each face
% pg         A series of six axis handles to the grids for each face
% pc         A series of six axis handles to the continents being plotted
%
% SEE ALSO: PLOTONCUBE, PLOTONCHUNK, READPETSCBINARYVEC
% 
% Last modified by fjsimons-at-alum.mit.edu, 05/20/2010

% Define defaults
defval('vname','vecx.petsc')
defval('dshell',27)
defval('tipe','2D')
defval('kind',1)
defval('vscale',[-1 1])
defval('colmap','kelicol')
defval('phil',1)
defval('setnan',1)

% This here needs to be had from cubeunfolded.pdf
load(fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS','radii_129'))
% Now downsample this to 37 layers
defval('depths',radii(129)-...
       [radii([2 4 7:4:95 98 100 103 107 110 112 115:4:123 126])...
	radii(127)+[radii(128)-radii(127)]/2 ...
	radii(128)+[radii(129)-radii(128)]/2])

if isstr(vname)
  % The input is read by Sergey Voronin's mex code READPETSCBINARYVEC
  mexpath=fullfile('/u/fjsimons/STUDENTS/SergeyVoronin',...
		   'gaussian_balls_analysis_fjs/mex_files/');
  addpath(mexpath);
  datapath=fullfile('/u/fjsimons/STUDENTS/SergeyVoronin',...
		    'gaussian_balls_analysis_fjs/data');
  
  if exist('readPetscBinaryVec.mexa64','file')~=3
    % The mex file is being recompiled
    disp(sprintf('Compiling %s',fullfile(mexpath,'readPetscBinaryVec.c')))
    mex(fullfile(mexpath,'readPetscBinaryVec.c'))
    disp('Make sure the mex file ends up somewhere reasonable')
  end
  
  % Read the solution in linear fashion
  x=readPetscBinaryVec(fullfile(datapath,vname));
  % Write it? writePetscBinaryVec('x',fullfile(datapath,'bla'));
  
  % Now rearrange over the depths and the chunks
  nchunk=6;
  ndepth=37;
  NND=length(x)/6;
  [Nxi,Neta]=deal(sqrt(NND/ndepth));
  
  % Initialize the vector at the shell that will be plotted
  v=zeros(Nxi,Neta,nchunk);
  
  % Pick out the chunks at the right depth
  for index=1:nchunk
    %   v(:,:,index)=reshape(x([Nxi*Neta*(dshell-1)+1:Nxi*Neta*dshell]...
    % 			 +(index-1)*NND)',...
    % 		       Nxi,Neta);
    v(:,:,index)=reshape(x((dshell:ndepth:NND)+(index-1)*NND),Nxi,Neta);
  end
elseif isstruct(vname)
  fnX=fieldnames(vname);
  szX=size(vname.(fnX{1}));
  nchunk=length(fnX);
  ndepth=szX(3);
  
  % Initialize the vector at the shell that will be plotted
  v=zeros(szX(1),szX(2),nchunk);

  % Rearrange the structure so that it is well presented
  for index=1:nchunk    
    v(:,:,index)=tindeks(vname.(fnX{index}),dshell);
  end
end

% Check we're dealing with the right depths here
difer(ndepth-length(depths),[],0,NaN)

clf
ah=gca;
% Now I can feed this straight to PLOTONCUBE
[p,pg]=plotoncube(v,tipe,kind,[],[],[],vscale,colmap,phil,setnan);

if strcmp(tipe,'2D')
  % Plot the continents on top
  [a,pc]=plotcont([],[],9);
  % Add a nice colorbar
  colpos=[0.5616    0.1714    0.3143    0.0298];
  if kind==1
    [cb,xcb]=addcb(colpos,vscale,vscale,colmap,range(vscale)/5,0);
    set(cb,'fonts',8)
    set(xcb,'string','anomaly','fonts',8)
  else kind==2
    cb=colorbarf('hor',8,'Helvetica',colpos);
    set(gcf,'NextPlot','add')
    axes(cb); xlabel('anomaly','fonts',8)
  end
  axes(ah)
end

% Add a good legend
if isstr(vname)
  kname=suf(datapath,'/');
  t(1)=title(sprintf('%s/%s',kname,nounder(vname)));
elseif isstruct(vname)
  % Change the vname to the inputname in the hereafter
  vname=inputname(1);
  t(1)=title(sprintf('The structure variable %s',vname));
end
xloc=-0.75;
t(2)=text(xloc,3.50,...
    sprintf('layer %i / % i',dshell,ndepth));
t(3)=text(xloc,3.25,...
    sprintf('depth %8.3f km',depths(dshell)));
t(4)=text(xloc,3.00,...
	  sprintf('min = %6.3f ; max = %6.3f',...
	    min(v(:)),max(v(:))));
% Read number of nonzeros in the wavelet domain
defval('nnz',NaN)
t(5)=text(xloc,2.75,...
	  sprintf('nnzw = %6.3i',nnz));
defval('wavlet','db4')
t(6)=text(xloc,2.50,...
	  sprintf('wavelet = %s',wavlet));
defval('algo','ISD/FISTA')
t(6)=text(xloc,2.25,...
	  sprintf('algorithm = %s',algo));
defval('edg','on')
t(7)=text(xloc,2.00,...
	  sprintf('edge = %s',edg));
defval('precon','on')
t(8)=text(xloc,1.75,...
	  sprintf('precon = %s',precon));
defval('mra','on')
t(9)=text(xloc,1.50,...
	  sprintf('mra = %s',mra));
defval('dfit',NaN)
t(10)=text(xloc,1.25,...
	  sprintf('datafit = %s',dfit));
defval('scalingthreshold','off')
t(11)=text(xloc,1.00...
	  sprintf('thresh scaling = %s',scalingthreshold));

keyboard

% Cosmetics
fig2print(gcf,'portrait')
if any(abs(vname)==46)
  vname=pref(vname);
end
figdisp([],sprintf('%s_%2.2i_%i',vname,dshell,ndepth),[],0)

% Prepare output
vals={p,pg,pc};
varargout=vals(1:nargout);
