function viewsg(ds, i, fignum)
% VIEWSG(ds, i, fignum) - View the ith synaptogram in Figure fignums

if(~isfield(ds,'sg'))
    sg = loadsg(ds, i);
else
    sg = ds.sg(i);
end

if(~exist('fignum', 'var'))
    fignum = gcf;
end

% title
t =  sprintf('Synaptogram for Synapse %d [ %s ]', i, sg.coordstr);

vspacing = 1;
im = [];

for c = 1:ds.nch
    row = viewsgrow(sg.im(c,:,:,:), [0 0 0], vspacing);
    sz = size(row);
    im = [im; row];
    
end

figure(fignum), clf;
set(gcf, 'Position', [827   476   817   669]);
set(gcf, 'Name', t);
set(gcf, 'NumberTitle', 'off');
set(gcf, 'Color', [0 0 0]);

imagesc(max(im,0));
yticks = vspacing + ds.sgdim(2)/2 + 1 : ds.sgdim(2) + 2*vspacing : sz(2);
set(gca, 'YColor', 'k');
axorig = gca;
set(gca, 'XTick', []);

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
   'XColor','k','YColor','w', 'XTick', []);

% left channel labels
set(ax1, 'YLim', get(axorig, 'YLim'));      
set(ax1, 'YTick', yticks);
set(ax1, 'YTickLabel', fliplr(ds.chlist));

% right channel labels
set(ax2, 'YLim', get(axorig, 'YLim'));
set(ax2, 'YTick', yticks);
set(ax2, 'YTickLabel', fliplr(ds.chlist));

sectlabels = -floor(ds.sgdim(3)/2):floor(ds.sgdim(3)/2);

set(ax1, 'XLim', get(axorig, 'XLim'));
set(ax1, 'XTick', yticks);
set(ax1, 'XTickLabel', sectlabels);
xlabel(ax1,'Section Plane')

title(ax1, t, 'Color', 'w');

% black overlay
axov1 = axes('Position',get(gca,'Position'),...
   'XAxisLocation','bottom',...
   'YAxisLocation','left',...
   'Color','none',...
   'XColor','k','YColor','k', ...
   'XLim', get(ax1, 'XLim'), 'XTick', get(ax1,'XTick'), 'XTickLabel', {}, 'TickDir', 'out', ...
   'YTick', yticks, 'YLim', get(axorig, 'YLim'), 'YTickLabel', {});

axov2 = axes('Position',get(gca,'Position'),...
   'XAxisLocation','bottom',...
   'YAxisLocation','right',...
   'Color','none',...
   'XColor','k','YColor','k', 'XTick', [], 'TickDir', 'in', ...
   'YTick', yticks, 'YLim', get(axorig, 'YLim'), 'YTickLabel', {});





