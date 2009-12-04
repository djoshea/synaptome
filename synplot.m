function synplot(ds, features, errors)

if(~exist('features', 'var'))
    features = {'IB_VGlut1', 'IB_PSD95'};
end
if(~exist('errors', 'var'))
    errors = zeros(ds.ntrain,1);
end

t = sprintf('%s vs. %s', features{1}, features{2});
[ft idx] = getfeature(ds, features);
    
clf; hold on
% set(gcf, 'Position', [244   543   560   420]);
set(gcf, 'Name', t);
set(gcf, 'NumberTitle', 'off');
set(gcf, 'defaulttextinterpreter', 'none');

h = zeros(size(ds.labelnames));
for ci = 1:length(ds.labelnames)
    h(ci) = line(0, 0, 'Marker', 'o', 'MarkerSize', 8, 'LineStyle', 'none', ...
        'MarkerFaceColor', ds.labelcolors(ci,:), ...
        'MarkerEdgeColor', ds.labelcolors(ci,:));
end
legend(ds.labelnames, 'Location', 'Best');
delete(h);

for i = 1:ds.ntrain
    if(~ds.trainactive(i))
        continue;
    end
       
    ci = ds.trainlabel(i);
    
    % mark synapse as error?
    if(errors(i))
        h = line(ft(i,1), ft(i,2));
        set(h, 'Marker','x', 'LineStyle', 'none', ...
        'Color', ds.labelcolors(ci,:), ...
        'MarkerFaceColor', ds.labelcolors(ci,:), ...
        'MarkerEdgeColor', ds.labelcolors(ci,:), ...
        'MarkerSize', 15, ...
        'Tag', num2str(i));
    end
    
    h = line(ft(i, 1), ft(i, 2));
    set(h, 'Marker','o', 'LineStyle', 'none', ...
     'Color', ds.labelcolors(ci,:), ...
     'MarkerFaceColor', ds.labelcolors(ci,:), ...
     'MarkerEdgeColor', ds.labelcolors(ci,:), ...
     'MarkerSize', 4, ...
     'ButtonDownFcn', @callback, ...
     'Tag', num2str(i));
 
 
        
end
set(gcf, 'UserData', '');

xlabel(features{1});
ylabel(features{2});
title(t);

function callback(src, event)
   
tag = str2num(get(src,'Tag'));
x = get(src, 'XData');
y = get(src, 'YData');

hmark = get(gcf, 'UserData');
if(~strcmp(hmark, ''))
    delete(hmark);
end

hcurrent = plot(x,y, 'Marker', 'o', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'g');
set(gcf,'UserData', hcurrent);

viewsg(ds,tag, 99);

end

end