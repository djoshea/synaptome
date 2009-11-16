%% Profile: All Synapses

ds.labelnames = {'none', 'Glut', 'GABA', 'both'};

figure(1), clf;
set(gcf, 'NumberTitle', 'off')
set(gcf, 'Name', 'Synaptic Density Profile: All Synapses');
[ax h1 h2] = plotyy(ds.dbins(1:end-2)*ds.dxy, ds.dcounts(1,1:end-2) / ds.dbinvol, ...
                    ds.dbins(1:end-2)*ds.dxy, ds.dcounts(2,1:end-2) / ds.dbinvol);
hold off
hleg = legend(ds.labelnames{ds.showlabels}, 'Location', 'Best');
set(hleg, 'FontSize', 12);
set(hleg, 'Position', [0.6176    0.2093    0.0967    0.2682]);
box off
legendboxoff
xlabel('Cortical Depth from Pia ($\mu m$)');
set(ax(1), 'XLim', [0 ds.dbins(end-2)*ds.dxy]);
set(ax(2), 'XLim', [0 ds.dbins(end-2)*ds.dxy]);
set(ax(1), 'YTick', 0:0.4:2);
set(ax(2), 'YTick', 0:0.02:0.1);
set(ax(1), 'YColor', [1 0 0]);
set(ax(2), 'YColor', [0 0 1]);
set(get(ax(1), 'YLabel'), 'String', 'Glut Synapses per $\mu m^3$');
set(get(ax(2), 'YLabel'), 'String', 'GABA Synapses per $\mu m^3$');
set(h1, 'Color', [1 0 0]);
set(h2, 'Color', [0 0 1]);
title('Cortical Synaptic Density Profile');
set(gcf, 'Position', [100 10 1050 185]);

% %% Profile YFP and Synapses onto YFP
% figure(2), clf;
% set(gcf, 'NumberTitle', 'off')
% set(gcf, 'Name', 'Synaptic Density Profile: YFP');
% [ax h1 h2] = plotyy(ds.dbins(1:end-2)*ds.dxy, ds.mask{2}.dcounts(1,1:end-2) / ds.dbinvol, ...
%                     ds.dbins(1:end-2)*ds.dxy, ds.mask{2}.dcounts(2,1:end-2) / ds.dbinvol);
% hold on
% 
% yfpdens = ds.mask{1}.dpxcounts(1:end-2) / ds.dbinvol;
% yfpnorm = yfpdens * max(get(ax(1),'YLim')) / max(yfpdens);
% plot(ds.dbins(1:end-2)*ds.dxy, yfpnorm, 'g-', 'LineWidth', 2);
% 
% legend({'Glut', 'GABA', 'YFP'}, 'Location', 'Best');
% hold off
% box off
% legendboxoff
% xlabel('Cortical Depth from Pia ($\mu m$)');
% set(ax(1), 'XLim', [0 ds.dbins(end-2)*ds.dxy]);
% set(ax(2), 'XLim', [0 ds.dbins(end-2)*ds.dxy]);
% set(ax(1), 'YTick', 0:0.4:2);
% set(ax(2), 'YTick', 0:0.02:0.1);
% set(ax(1), 'YColor', [1 0 0]);
% set(ax(2), 'YColor', [0 0 1]);
% set(get(ax(1), 'YLabel'), 'String', 'Synapses per $\mu m^3$');
% set(get(ax(2), 'YLabel'), 'String', 'Synapses per $\mu m^3$');
% set(h1, 'Color', [1 0 0]);
% set(h2, 'Color', [0 0 1]);
% title('YFP Synaptic Density Profile');

%% Profile YFP, Synapses onto YFP, Peters' Rule Predicted Syn onto YFP

yfpdensmask = 2;
synmask = 2;

% pixels / um^3
yfpdens = ds.mask{yfpdensmask}.dpxcounts(1:end-2) * (ds.dxy^2 * ds.dz) / ds.dbinvol;
exp_maskdens = repmat(yfpdens,2,1) .* ds.dcounts(:,1:end-2) / ds.dbinvol;

% pixels / total pixels
yfpprop = double(yfpdens) / 20; 

figure(3), clf;
set(gcf, 'NumberTitle', 'off')
set(gcf, 'Name', 'Synaptic Density Profile: YFP');
[ax h1 h2] = plotyy(ds.dbins(1:end-2)*ds.dxy, ds.mask{synmask}.dcounts(1,1:end-2) / ds.dbinvol, ...
                    ds.dbins(1:end-2)*ds.dxy, ds.mask{synmask}.dcounts(2,1:end-2) / ds.dbinvol);
h1e = line(ds.dbins(1:end-2)*ds.dxy, exp_maskdens(1,:), 'Parent', ax(1)); 
h2e = line(ds.dbins(1:end-2)*ds.dxy, exp_maskdens(2,:), 'Parent', ax(2));

set(ax(1), 'YLim', [0 0.5]);
set(ax(2), 'YLim', [0 0.025]);

yfpdens = ds.mask{yfpdensmask}.dpxcounts(1:end-2) / ds.dbinvol;
yfpnorm = yfpdens * max(get(ax(1), 'YLim')) / max(yfpdens);

hy = line(ds.dbins(1:end-2)*ds.dxy, yfpdens, 'Parent', ax(1));
set(hy, 'Color', [0 1 0]);
set(hy, 'LineStyle', '-');
set(hy, 'LineWidth', 2);

hleg = legend({'Glut', 'Glut Expected', 'YFP Density (A.U.)', 'GABA', 'GABA Expected'}, 'Location', 'Best');
set(hleg, 'FontSize', 12);
set(hleg, 'Position', [0.1760    0.3317    0.1652    0.5811]);

hold off
box off
legendboxoff
xlabel('Cortical Depth from Pia ($\mu m$)');
set(ax(1), 'XLim', [0 ds.dbins(end-2)*ds.dxy]);
set(ax(2), 'XLim', [0 ds.dbins(end-2)*ds.dxy]);
set(ax(1), 'YTick', 0:0.1:0.5);
set(ax(2), 'YTick', 0:0.005:0.025);
set(ax(1), 'YColor', [1 0 0]);
set(h1e, 'Color', [1 0 0]);
set(h1e, 'LineWidth', 1);
set(h1e, 'LineStyle', '--');
set(h2e, 'Color', [0 0 1]);
set(h2e, 'LineWidth', 1);
set(h2e, 'LineStyle', '--');

set(ax(2), 'YColor', [0 0 1]);
set(get(ax(1), 'YLabel'), 'String', 'Glut Synapses per $\mu m^3$');
set(get(ax(2), 'YLabel'), 'String', 'GABA Synapses per $\mu m^3$');
set(h1, 'Color', [1 0 0]);
set(h2, 'Color', [0 0 1]);
title('YFP Synaptic Density Profile');

h1ov = line(ds.dbins(1:end-2)*ds.dxy, ds.mask{synmask}.dcounts(1,1:end-2) / ds.dbinvol, 'Parent', ax(1));
h2ov = line(ds.dbins(1:end-2)*ds.dxy, ds.mask{synmask}.dcounts(2,1:end-2) / ds.dbinvol, 'Parent', ax(2));
h1eov = line(ds.dbins(1:end-2)*ds.dxy, exp_maskdens(1,:), 'Parent', ax(1)); 
h2eov = line(ds.dbins(1:end-2)*ds.dxy, exp_maskdens(2,:), 'Parent', ax(2));

set(h1eov, 'Color', [1 0 0]);
set(h1eov, 'LineWidth', 1);
set(h1eov, 'LineStyle', '--');
set(h2eov, 'Color', [0 0 1]);
set(h2eov, 'LineWidth', 1);
set(h2eov, 'LineStyle', '--');

set(h1ov, 'Color', [1 0 0]);
set(h2ov, 'Color', [0 0 1]);

set(gcf, 'Position', [100 500 1250 185]);
