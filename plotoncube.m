function varargout=...
    plotoncube(v,tipe,kind,x,y,z,vscale,colmap,phil,setnan,wgrid)
% [h,pg,txtywk,pgw]=...
%   PLOTONCUBE(v,tipe,kind,x,y,z,vscale,colmap,phil,setnan,wgrid)
%
% Plots SINGLE-DEPTH cubed-sphere data.
%
% What needs to be checked is that we used to have an odd number of
% points, redundant at the edges, and now we have an even number of
% points, which presumably define pixel centers instead of nodes.
%
% INPUT:
%
% v       The NxNx6 values to plot on the cubed sphere at a single depth.
%         The first dimension v(i,:,:) corresponds to xi and points up in
%         a two-dimensional rendition; it is the original longitude of
%         the x+ chunk. The second dimension v(:,i,:) corresponds to eta
%         and points across in a two-dimensional plot; it is the original
%         colatitude of the x+ chunk. The third dimension is the face
%         number. For the definition of the faces, see CUBEMATS.
% tipe    '3D' [default] for a 3-D spherical rendering
%         '2D' for a 2-D Cartesian rendering
% kind    1 for a 3-D plot using SURF on a globe [default]
%         2 for a 3-D plot using a relief SURFACE on a globe
%         1 for a 2-D plot using IMAGEFNAN [default]
%         2 for a 2-D plot using IMAGESC
%         3 for a 2-D plot using SURF
% x,y,z   The NxNx6 coordinates of the cubed sphere [defaulted]
% vscale  The color saturation in case you want this fixed
% colmap  A string name with the color map [default: 'kelicol']
% phil    Flag to reorder the faces (1) or not (0)
% setnan  1 Set values smaller than 100*eps to NaN [default]
%         0 Don't 
%         A scalar that defines the threshold, but see also SETNANS
% wgrid   In two-dimensional mode, also uses FRIDPLOTW for the grid using
%         wgrid.N and wgrid.J if they are defined
%
% OUTPUT:
%
% h       A series of six axis handles to the plots for each face
% pg      A series of six axis handles to the grids for each face
% txtyw   A series of six x,y,width coordinates in the top left
% pgw     A series of six axis handles to the wavelet grids for each face
%
% SEE ALSO: PLOTONCHUNK, PLOTONCUBE2, PLOTONCUBE3, PLOTCONT
% 
% EXAMPLE:
%
% plotoncube('demo1') % With letters in 2D (screen-resolution dependent!)
% plotoncube('demo2') % With letters in 3D (screen-resolution dependent!)
% plotoncube('demo3') % With blank faces in 2D
% plotoncube('demo4') % With blank faces in 3D
%
% Last modified by fjsimons-at-alum.mit.edu, 2/16/2011

if ~isstr(v)
  % Define defaults
  defval('tipe','3D')
  defval('kind',1)
  defval('phil',0)
  defval('colmap','kelicol')
  defval('setnan',1)
  defval('wgrid',[])

  if phil
    % Change to the new definition of the faces
    v=v(:,:,[5 4 6 2 3 1]);
    disp('Input was in Judd/Vetter chunk numbering scheme')
  end
  
  % Construct the sphere if you haven't already
  isod=mod(length(v(:,:,1)),2);
  lfin=nextpow2(size(v(:,:,1),1))-isod;
  % The double logical operators are crucial here in this case
  if strcmp(tipe,'3D') && [~exist('x') || [exist('x') & isempty(x)]]
    disp('This needs to be adapted to only return 128 points')
    [x,y,z]=cube2sphere(lfin);
  end

  % Determine a common size and caxis for all faces
  N=2^lfin+1;
  defval('vscale',halverange(v,75,NaN));

  % Now do the actual plotting - also allow for a single pane
  for in=1:min(6,size(v,3))
    switch tipe
     case '3D'
      switch kind
       case 1
	v(v<vscale(1))=vscale(1);
	v(v>vscale(2))=vscale(2);
	p{in}=surf(x(:,:,in),y(:,:,in),z(:,:,in),flipud(v(:,:,in)));
       case 2
	% "Vertical" exaggeration
	defval('rang',0.05)
	% Create the surface
	p{in}=surface(x(:,:,in),y(:,:,in),z(:,:,in),...
		      'FaceColor','texture','Cdata',flipud(v(:,:,in)));
	% Rescale the data, protect against NaN
	vs=flipud(v(:,:,in));
	% Keep the zero-level explicitly in there - don't do scale
	vs=1+rang*vs;
	xs=x(:,:,in).*vs;      
	ys=y(:,:,in).*vs;
	zs=z(:,:,in).*vs;
	set(p{in}, 'xdata',xs,'ydata',ys,'zdata',zs)
	view(140,30) % So x and y face at you
	set(gca, 'dataaspectratio', [1 1 1], 'cameraviewangle', 7)
      end
      hold on
      pg=NaN;
      txtyw=NaN;
     case '2D'
      % The plot midpoints referenced to the angles/size of the cube
      xup=floor(in/2)*pi/2;      % xup=floor(in/2)*N;
      yup=(ceil(in/2)-1)*pi/2;   % yup=(ceil(in/2)-1)*N;
      % Identify the coordinates of the top left and bottom right
      c11=[xup-pi/4 yup+pi/4];   % c11=[xup+1 yup+1];
      cmn=[xup+pi/4 yup-pi/4];   % cmn=[xup+N yup+N];
      % Slight variation for the edges in each coordinate direction
      xups=xup+[-pi/4  pi/4];    % xups=xup+[1 N];
      yups=yup+[ pi/4 -pi/4];    % yups=yup+[1 N];
      % What is being plotted? First, flip upside down
      vface=flipud(v(:,:,in));

      % Save for later use
      txtyw(in,:)=[xups(1) yups(1) range(xups) range(yups)];

      switch kind
       case 1
	p{in}=imagefnan(c11,cmn,vface,colmap,vscale,[],0,setnan); hold on
	pg{in}=fridplot(xups,yups);
	axis xy image
        if ~isempty(wgrid)
          pgw{in}=fridplotw(wgrid.N,wgrid.J,[xups yups]); 
        end
       case 2
	p{in}=imagesc(xups,yups,vface); hold on
	pg{in}=fridplot(xups,yups);
	colormap(colmap)
	caxis(vscale)
	axis xy image
        if ~isempty(wgrid)
          pgw{in}=fridplotw(wgrid.N,wgrid.J,[xups yups]); 
        end
       case 3
	p{in}=surf(linspace(xups(1),xups(2),N),...
		   linspace(yups(1),yups(2),N),...
 		   vface); hold on
	colormap(colmap)
	shading flat
	view(10,80)
	pg=NaN;
      end
    end
  end
  hold off
  axis off

  defval('pgw',[])
  % Prepare output
  vals={p,pg,txtyw,pgw};
  varargout=vals(1:nargout);
elseif strcmp(v,'demo1')
  % Two-dimensional plots
  v=nan(65,65,6); 
  for in=1:6; v(:,:,in)=flipud(imageletter(in)); v(1:5,1:10,in)=1; end
  clf
  [p,pg]=plotoncube(v,'2D',2);
  figdisp([],1,[],1)
elseif strcmp(v,'demo2')
  % Three-dimensional plots
  v=nan(65,65,6); 
  for in=1:6; v(:,:,in)=flipud(imageletter(in)); v(1:5,1:10,in)=1; end
  clf
  [p,pg]=plotoncube(v,'3D'); shading faceted; colormap(flipud(gray))
  axis on
  xlabel('x'); ylabel('y'); zlabel('z')
  view(64,32)
  hold on
  plotcont([],[],3)
elseif strcmp(v,'demo3')
  % Two-dimensional plots
  v=nan(65,65,6); for in=1:6; v(1:5,1:10,in)=10; end
  clf
  [p,pg]=plotoncube(v,'2D');
  figdisp([],3,[],1)
elseif strcmp(v,'demo4')
  % Three-dimensional plots
  v=nan(65,65,6); 
  for in=1:6; 
    v(1:5,1:5,in)=10; 
    v(1:3,:,in)=1; v(end-2:end,:,in)=1; 
    v(:,1:3,in)=1; v(:,end-2:end,in)=1;     
  end
  clf
  [p,pg]=plotoncube(v,'3D');
  figdisp([],4,[],1)
  shading flat
end
