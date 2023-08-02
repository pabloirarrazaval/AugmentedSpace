function pretty_plot(varargin)
% pretty_plot([data]). Makes prettier the current axes.
%
% data is an optional argument with all or some of these fields:
% data.axis_color = [r g b] (defaults to [0.5 0.5 0.5])
% data.axis_width = w (defaults to 1)
% data.axis_centered = false|true (defautls to true: puts axis as close
%                           as posible to the origen, otherwise put them to
%                           the left or bottom)
% data.xlabel = string (will be interpreted as latex if wrapped in $)
% data.xticks = array of string's cells (add tick marks at thos positions)
% data.ylabel = string (will be interpreted as latex if wrapped in $)
% data.yticks = array of string's cells (add tick marks at thos positions)
% data.fontsize = n (defaults to 14)
% data.font_color = [r g b] (defaults to [0.3 0.3 0.3])
%
% Pablo Irarrazaval (July 2023)

% Defaults
axis_color = [0.5 0.5 0.5];
axis_width = 1;
axis_centered = true;
xlabel = [];
ylabel = [];
xticks = {};
yticks = {};
fontsize = 14;
font_color =[0.3 0.3 0.3];
if nargin > 0 % options are specfied (extra inputs are ignored)
    option = varargin{1};
    if ~isstruct(option)
        error('Input must be a structure with data.');
    end
    if isfield(option,'axis_color'), axis_color = option.axis_color; end
    if isfield(option,'axis_width'), axis_width = option.axis_width; end
    if isfield(option,'axis_centered'), axis_centered = option.axis_centered; end
    if isfield(option,'xlabel'), xlabel = option.xlabel; end
    if isfield(option,'xticks'), xticks = option.xticks; end
    if isfield(option,'ylabel'), ylabel = option.ylabel; end
    if isfield(option,'yticks'), yticks = option.yticks; end
    if isfield(option,'fontsize'), fontsize = option.fontsize; end
    if isfield(option,'font_color'), font_color = option.font_color; end
end

% Expected behavior compatible with hold on/off (keeps data)
% Clear annotations (axes in our case) and labels
delete(findall(gcf,'type','annotation'))
delete(findall(gcf,'Type','text'))

% White background
set(gcf,'Color','white');

% Select the current axes
ax = gca;

% Erase current axes
axis('off'); % Erases all
% ax.Box = 'off'; % Erases the box
% ax.XTick = []; % Erases ticks
% ax.YTick = [];

lim = axis; % actual data limits

% We will use normalized units
originalunits = normalized_units(ax);

axpos = ax.Position; % normalized position of axes

%% Draw axes
% Compute position for drawing axes
limwidth = lim(2)-lim(1);
limheight = lim(4)-lim(3);
if axis_centered
    cxd = min(max(0,lim(1)),lim(2)); % x position in data units
    cyd = min(max(0,lim(3)),lim(4)); % y position in data units
    cx = axpos(1) + (cxd - lim(1))/limwidth*axpos(3); % center in x
    cy = axpos(2) + (cyd - lim(3))/limheight*axpos(4); % center in y
else
    cxd = lim(1); % x position in data units
    cyd = lim(3); % y position in data units
    cx = axpos(1);
    cy = axpos(2);
end

% Draw arrows in axes positions
% Extra length so that the arrow is further away
xend = 1-(1-axpos(1)-axpos(3))/2;
anx = annotation('arrow',[axpos(1) xend],[cy cy]); % x axis
yend = 1-(1-axpos(2)-axpos(4))/2;
any = annotation('arrow',[cx cx],[axpos(2) yend]); % y axis
% Changes width and color
set(anx,'LineWidth',axis_width,'Color',axis_color);
set(any,'LineWidth',axis_width,'Color',axis_color);

%% Write labels in normalized units
% Normalized units for text are scaled to InnerPosition of axes
ctx = (cx-axpos(1))/axpos(3); % center of axes in text coordinates
cty = (cy-axpos(2))/axpos(4);
if ~isempty(xlabel)
    if xlabel(1)=='$', interp = 'latex'; else, interp = 'none'; end
    x = 1.08; % a small percentage to the right
    y = cty - 0.015; % a small percentage down
    t = text(x,y,xlabel,'Interpreter',interp,'Units','normalized',...
        'VerticalAlignment','top','HorizontalAlignment','right');
    set(t,'FontSize',fontsize,'Color',font_color);
end
if ~isempty(ylabel)
    if ylabel(1)=='$', interp = 'latex'; else, interp = 'none'; end
    x = ctx - 0.015; % a small percentage to the left
    y = 1.08; % a small percentage up
    t = text(x,y,ylabel,'Interpreter',interp,'Units','normalized',...
        'VerticalAlignment','top','HorizontalAlignment','right');
    set(t,'FontSize',fontsize,'Color',font_color);
end

%% Write ticks
Nticks = length(yticks); % in Y
for n = 1:Nticks
    y = str2num(yticks{n});
    x = cxd-limwidth*0.01;
    t = text(x,y,yticks{n},...
        'VerticalAlignment','middle','HorizontalAlignment','right');
    set(t,'FontSize',12,'Color',font_color);
end
Nticks = length(xticks); % in X
for n = 1:Nticks
    x = str2num(xticks{n});
    y = cyd;
    t = text(x,y,xticks{n},...
        'VerticalAlignment','top','HorizontalAlignment','center');
    set(t,'FontSize',12,'Color',font_color);
end

%%
restore_units(ax);

%% Auxiliary functions
    function originalunits = normalized_units(ax)
        % Changes graph to Normalized units to be sure all dimensions work
        originalunits = get(ax, 'Units'); % save current units to restore them later
        set(ax, 'Units', 'Normalized');
    end

    function restore_units(ax)
        % Go back to old units
        set(ax, 'Units', originalunits);
    end

end

