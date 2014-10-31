function handle=colorbarf(varargin)
% colorbarf.m
% Displays colorbar axis.
% Takes four input arguments:  LOCATION ('hor', 'ver')
%                              FONTSIZE
%                              FONTNAME ('Times, 'Helvetica', and so on)
%                 (optional)   POSITION ([Left Bottom Width Height])
% If using POSITION, must also specify FONTSIZE, FONTNAME
% From the original colorbar.m, modified by FJS, June 27, 1997.
% AFTER THIS MAY WANT TO DO set(gcf,'NextPlot','add')

%COLORBAR Display color bar (color scale).
%   COLORBAR('vert') appends a vertical color scale to the current
%   axis. COLORBAR('horiz') appends a horizontal color scale.
%
%   COLORBAR(H) places the colorbar in the axes H. The colorbar will
%   be horizontal if the axes H width > height (in pixels).
%
%   COLORBAR without arguments either adds a new vertical color scale
%   or updates an existing colorbar.
%
%   H = COLORBAR(...) returns a handle to the colorbar axis.

%   Clay M. Thompson 10-9-92
%   Copyright (c) 1984-96 by The MathWorks, Inc.
%   $Revision: 5.18 $  $Date: 1996/10/22 15:10:46 $

%   If called with COLORBAR(H) or for an existing colorbar, don't change
%   the NextPlot property.
changeNextPlot = 1;

loc=varargin{1};
ftsize=varargin{2};
ftname=varargin{3};

if nargin==4 
pospos=varargin{4};
end

if nargin<1, loc = 'vert'; end
ax = [];
if nargin==1,
    if ishandle(loc)
        ax = loc;
        if ~strcmp(get(ax,'type'),'axes'),
            error('Requires axes handle.');
        end
        units = get(ax,'units'); set(ax,'units','pixels');
        rect = get(ax,'position'); set(ax,'units',units)
        if rect(3) > rect(4), loc = 'horiz'; else loc = 'vert'; end
        changeNextPlot = 0;
    end
end

% Determine color limits by context.  If any axes child is an image
% use scale based on size of colormap, otherwise use current CAXIS.

ch = get(gca,'children');
hasimage = 0; t = [];
cdatamapping = 'direct';

for i=1:length(ch),
    typ = get(ch(i),'type');
    if strcmp(typ,'image'),
        hasimage = 1;
        cdatamapping = get(ch(i), 'CDataMapping');
    elseif strcmp(typ,'surface') & ...
            strcmp(get(ch(i),'FaceColor'),'texturemap') % Texturemapped surf
        hasimage = 2;
        cdatamapping = get(ch(i), 'CDataMapping');
    elseif strcmp(typ,'patch') | strcmp(typ,'surface')
        cdatamapping = get(ch(i), 'CDataMapping');
    end
end

if ( strcmp(cdatamapping, 'scaled') )
    if hasimage,
        if isempty(t); 
            t = caxis; 
        end
    else
        t = caxis;
        d = (t(2) - t(1))/size(colormap,1);
        t = [t(1)+d/2  t(2)-d/2];
    end
else
    if hasimage,
        t = [1, size(colormap,1)]; 
    else
        t = [1.5  size(colormap,1)+.5];
    end
end

h = gca;

if nargin==0,
    % Search for existing colorbar
    ch = get(gcf,'children'); ax = [];
    for i=1:length(ch),
        d = get(ch(i),'userdata');
        if prod(size(d))==1 & isequal(d,h), 
            ax = ch(i); 
            pos = get(ch(i),'Position');
            if pos(3)<pos(4), loc = 'vert'; else loc = 'horiz'; end
            changeNextPlot = 0;
            break; 
        end
    end
end

origNextPlot = get(gcf,'NextPlot');
if strcmp(origNextPlot,'replacechildren') | strcmp(origNextPlot,'replace'),
    set(gcf,'NextPlot','add')
end

if loc(1)=='v', % Append vertical scale to right of current plot
    
    if isempty(ax),
        units = get(h,'units'); set(h,'units','normalized')
        pos = get(h,'Position'); 
        [az,el] = view;
        stripe = 0.075; edge = 0.02; 
        if all([az,el]==[0 90]), space = 0.05; else space = .1; end
           set(h,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
        if exist('pospos')==1
           set(h,'Position',pos)
           pos=pospos;
	end      

        rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
        
        % Create axes for stripe
        ax = axes('Position', rect);
        set(h,'units',units)
    else
        axes(ax);
    end
    
    % Create color stripe
    n = size(colormap,1);
    image([0 1],t,(1:n)','Tag','TMW_COLORBAR'); set(ax,'Ydir','normal')
    set(ax,'YAxisLocation','right')
    set(ax,'xtick',[])
    
else, % Append horizontal scale to top of current plot
    
    if isempty(ax),
        units = get(h,'units'); set(h,'units','normalized')
        pos = get(h,'Position');
        stripe = 0.075; space = 0.1;
        set(h,'Position',...
            [pos(1) pos(2)+(stripe+space)*pos(4) pos(3) (1-stripe-space)*pos(4)])
        if exist('pospos')==1
           set(h,'Position',pos)
           pos=pospos;
	end      

        rect = [pos(1) pos(2) pos(3) stripe*pos(4)];
        
        % Create axes for stripe
        ax = axes('Position', rect);
        set(h,'units',units)
    else
        axes(ax);
    end
    
    % Create color stripe
    n = size(colormap,1);
    image(t,[0 1],(1:n),'Tag','TMW_COLORBAR'); set(ax,'Ydir','normal')
    set(ax,'ytick',[])    
end
set(ax,'userdata',h)
set(gcf,'CurrentAxes',h)

set(ax,'FontSize',[ftsize]')
set(ax,'FontName',ftname)

if changeNextPlot
    set(gcf,'Nextplot','ReplaceChildren')
end

if nargout>0, handle = ax; end
if exist('pospos')==1
set(ax,'Position',pospos)
end      
