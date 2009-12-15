function viewsg(ds, i, fignum)
% VIEWSG(ds, i, fignum) - View the ith synaptogram in Figure fignums

bgcol = [0.14 0.14 0.14];

if(~isfield(ds,'sg'))
    sg = loadsg(ds, i); % original 11^3 box of pixel values and some info
else
    sg = ds.sg(i);
end

if(~isfield(sg, 'str'))
    str = ds.labelnames{ds.trainlabel(i)};
else
    str = sg.str;
end

if(~exist('fignum', 'var'))
    fignum = 199;
end

% title
t =  sprintf('Synaptogram for Synapse %d [ %s ]', i, str);
vspacing = 1;
im = [];

for r = 1:length(ds.vis);
    rname = ds.vis{r};
    if(ischar(rname))
        dat = getchannel(ds,rname, i);
    else
        dat = [];
        dat.R = getchannel(ds,rname{1}, i);
        dat.G = getchannel(ds,rname{2}, i);
        dat.B = getchannel(ds,rname{3}, i);
    end
    row = viewsgrow(dat, bgcol, vspacing);
    sz = size(row);
    im = [im; row];
end

figure(fignum), clf;
% set(gcf, 'Position', [827   476   817   669]);
set(gcf, 'Name', t);
set(gcf, 'NumberTitle', 'off');
set(gcf, 'Color', bgcol);
set(gcf, 'ToolBar', 'none');

imagesc(max(im,0));
set(gca,'Units', 'normalized');
set(gca,'Position', [0.15 0.08 0.7 0.84]);
ystep = ds.sgdim(2) + 2*vspacing;
xstep = ds.sgdim(1) + 2*vspacing;
yticks = vspacing + ds.sgdim(2)/2 + 1 : ystep : length(ds.vis)*(ystep);
xticks = vspacing + ds.sgdim(1)/2 + 1 : xstep : vspacing + ds.sgdim(1)/2 +1 + xstep*ds.sgdim(3);
set(gca, 'YColor', bgcol);
axorig = gca;
set(gca, 'XTick', []);
set(gca, 'XColor', bgcol);

% left ylabels
ax1 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','left',...
           'Color','none',...
           'XColor','w','YColor','w', 'XTick', [], 'TickDir', 'out');
     
% right ylabels
ax2 = axes('Position',get(gca,'Position'),...
   'XAxisLocation','bottom',...
   'YAxisLocation','right',...
   'Color','none',...
   'XColor',bgcol,'YColor','w', 'XTick', []);

% left channel labels
set(ax1, 'YLim', get(axorig, 'YLim'));      
set(ax1, 'YTick', yticks);
set(ax1, 'YTickLabel', fliplr(ds.visname));

% right channel labels
set(ax2, 'YLim', get(axorig, 'YLim'));
set(ax2, 'YTick', yticks);
set(ax2, 'YTickLabel', fliplr(ds.visname));

sectlabels = -floor(ds.sgdim(3)/2):floor(ds.sgdim(3)/2);

set(ax1, 'XLim', get(axorig, 'XLim'));
set(ax1, 'XTick', xticks);
set(ax1, 'XTickLabel', sectlabels);
xlabel(ax1,'Section Plane')

title(ax1, t, 'Color', 'w');

% black overlay
axov1 = axes('Position',get(gca,'Position'),...
   'XAxisLocation','bottom',...
   'YAxisLocation','left',...
   'Color','none',...
   'XColor',bgcol,'YColor',bgcol, ...
   'XLim', get(ax1, 'XLim'), 'XTick', get(ax1,'XTick'), 'XTickLabel', {}, 'TickDir', 'out', ...
   'YTick', yticks, 'YLim', get(axorig, 'YLim'), 'YTickLabel', {});

axov2 = axes('Position',get(gca,'Position'),...
   'XAxisLocation','bottom',...
   'YAxisLocation','right',...
   'Color','none',...
   'XColor',bgcol,'YColor',bgcol, 'XTick', [], 'TickDir', 'in', ...
   'YTick', yticks, 'YLim', get(axorig, 'YLim'), 'YTickLabel', {});

syntable(ds,i);



