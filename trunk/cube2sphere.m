function varargout=cube2sphere(lfin,alfa,bita,gama,eo,sc)
% [x,y,z,J,N,dxi,deta,legs]=CUBE2SPHERE(lfin,alfa,bita,gama,eo,sc)
%
% Constructs a cubed sphere on 6 faces with 2^lfin[+1] nodes in each
% direction. Reference is Ronchi et al., doi:10.1006/jcph.1996.0047 
%
% INPUT:
%
% lfin      Number of subdivisions [default: 6]
% alfa      First Euler angle of wholesale tilt of all tiles [defaulted]
% beta      Second Euler angle of wholesale tilt of all tiles [defaulted]
% gama      Third Euler angle of wholesale tilt of all tiles [defaulted]
% eo        0 even number of points [default]
%           1 odd number of points
% sc        0 regular cubed sphere [default]
%           1 superchunk cubed sphere
%
% OUTPUT:
%
% x,y,z     The NxNx6 matrices with the Cartesian coordinates for the six faces
% J         The NxN Jacobian of the transformation for a single face
% N         The number of points in one grid dimension, i.e. 2^lfin+1
% dxi,deta  The spacing in xi and eta
% legs      A legend for the different chunk faces
%
% SEE ALSO: PLOTONCUBE, IMAGELETTER, SPHERE2CUBE, PLOTONCHUNK, CUBEMATS, CUBEJAC
% 
% EXAMPLE:
%
% cube2sphere('demo1')
% cube2sphere('demo2')
% cube2sphere('demo3')
%
% VOLUME CALCULATION:
% sum(6*dxi*deta*sum(sum(J(1:end-1,1:end-1))).*(rads/1000).^2*dr)
%
% Last modified by fjsimons-at-alum.mit.edu, 1/20/2011

% First define the defaults of all input variables
defval('lfin',6)
defval('alfa',[]);
defval('bita',[]);
defval('gama',[]);
defval('eo',0);
defval('sc',0);

if ~isstr(lfin)
  % Be in the know what the final dimension will be
  N=2^lfin+eo;

  % Get the rotation matrices
  [rottot,mats,legs]=cubemats(alfa,bita,gama);
  
  % Produce the Jacobian, the single-face coordinates and the angular spacing 
  [J,coordd,dxi,deta]=cubejac(N,N,sc);
  
  % Initialize the coordinate matrices
  [x,y,z]=deal(zeros(N,N,6));
  % Loop over all faces
  for f=1:6
    % Do the coordinate transform all at once, see SPHERE2CUBE for the inverse
    stuff=rottot*mats{f}*coordd;
        
    % Now distribute over the three three-dimensional vectors
    x(:,:,f)=reshape(stuff(1,:),N,N);
    y(:,:,f)=reshape(stuff(2,:),N,N);
    z(:,:,f)=reshape(stuff(3,:),N,N);
  end
  
  % Check that this is indeed still the unit sphere
  difer(x.^2+y.^2+z.^2-1,9,[],NaN)
  
  % The first and last rows/columns of J are exact repeats
  difer(J(1,:)-J(end,:),[],[],NaN)
  difer(J(:,1)-J(:,end),[],[],NaN)

  % Check that the integral over the surface is OK
  % The quality of this approximation depends on the number of
  % subdivisions carried out - see Ronchi eq. (20); deta*dxi is the
  % infinitesimal area right at the center of the block, where distortion
  % is minimal. But of course for finite subdivisions, this "rectangular"
  % area is not exactly to what it would be on the sphere
  %   disp(sprintf('Nominal area of sphere is %6.4f ; compare at %6.4f',...
  % 	       6*dxi*deta*sum(sum(J(1:end-1,1:end-1))),4*pi))
    
  % Provide output
  varns={x,y,z,J,N,dxi,deta,legs};
  varargout=varns(1:nargout);
elseif strcmp(lfin,'demo1')
  % Run the program
  [x,y,z,J]=cube2sphere; 
  % Define some colors
  cols={'c','m','y','k','r','b'};
  % Plot in succession
  for in=1:6
    p{in}=plot3(x(:,:,in),y(:,:,in),z(:,:,in),'o',...
		'MarkerF',cols{in},'MarkerE',cols{in});
    hold on; 
  end
  axis equal; hold off; axis off; set([p{:}],'MarkerS',3)
elseif strcmp(lfin,'demo2')
  % Run the program with no initial tilt
  [x,y,z,J,N]=cube2sphere(3,0,0,0); 
  % Plot in succession as a mesh
  for in=1:6
    p{in}=mesh(x(:,:,in),y(:,:,in),z(:,:,in),in*ones(N,N));
    hold on; 
  end
  colormap(gray(6))
  axis equal; hold off; axis off;
elseif strcmp(lfin,'demo3')
  clf
  % Run the program with no initial tilt
  [x,y,z,J,N,dx,de,legs]=cube2sphere(3,0,0,0); 
  % Plot if you like
  pols={'+','o','v','^','s','d'};
  for ind=1:6
    plot3(x(:,:,ind),y(:,:,ind),z(:,:,ind),...
	  'marker',pols{ind},'color','k')
    axis equal; axis([-1 1 -1 1 -1 1]);
    view(122,40); title(legs{ind})
    xlabel('x'); ylabel('y'); zlabel('z')
    pause
    hold on
  end
  hold off
end

