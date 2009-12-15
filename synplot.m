function synplot(ds, features, errors, twoclass)

if(~exist('features', 'var'))
    features = {'IB_VGlut1', 'IB_PSD95'};
end
if(~exist('errors', 'var'))
    errors = zeros(ds.ntrain,1);
else
    if(size(errors,1) ~= ds.ntrain)
        errors = zeros(ds.ntrain,1);
    end
end
if(~exist('twoclass', 'var'))
    dotwoclass = 0;
else
    dotwoclass = 1;
    twocolors = {'r', 'k'};
    Y = zeros(ds.ntrain,1);
    catnames = '';
    for i = 1:length(twoclass)
        idx = find(strcmp(ds.labelnames, twoclass{i}));
        if(numel(idx) == 0)
            error('Could not find class %s', twoclass{i});
        end
        if(i == 1)
            catnames = twoclass{i};
        else
            catnames = strcat(catnames, ', ', twoclass{i});
        end
        Y = Y | (ds.trainlabel == idx(1));
    end
end 
    
t = sprintf('%s vs. %s', features{1}, features{2});
[ft idx] = getfeature(ds, features);
    
clf; hold on
% set(gcf, 'Position', [244   543   560   420]);
set(gcf, 'Name', t);
set(gcf, 'NumberTitle', 'off');
set(gcf, 'defaulttextinterpreter', 'none');

if(~dotwoclass)
    for ci = 1:length(ds.labelnames)
        h = zeros(size(ds.labelnames));
        h(ci) = line(0, 0, 'Marker', 'o', 'MarkerSize', 8, 'LineStyle', 'none', ...
            'MarkerFaceColor', ds.labelcolors(ci,:), ...
            'MarkerEdgeColor', ds.labelcolors(ci,:));
    end
    legend(ds.labelnames, 'Location', 'Best');
else
    h = zeros(2,1);
    h(1) = line(0, 0, 'Marker', 'o', 'MarkerSize', 8, 'LineStyle', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    h(2) = line(0, 0, 'Marker', 'o', 'MarkerSize', 8, 'LineStyle', 'none', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    legend({catnames, 'Other'}, 'Location', 'Best');
end

delete(h);

for i = 1:ds.ntrain
    if(~ds.trainactive(i))
        continue;
    end
       
    if(dotwoclass)
        color = twocolors{2-Y(i)};
    else      
        ci = ds.trainlabel(i);
        color = ds.labelcolors(ci,:);
    end
    
    % mark synapse as error?
    if(errors(i))
        h = line(ft(i,1), ft(i,2));
        set(h, 'Marker','x', 'LineStyle', 'none', ...
        'Color', color, ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, ...
        'MarkerSize', 15, ...
        'Tag', num2str(i));
    end
    
    h = line(ft(i, 1), ft(i, 2));
    set(h, 'Marker','o', 'LineStyle', 'none', ...
     'Color', color, ...
     'MarkerFaceColor', color, ...
     'MarkerEdgeColor', color, ...
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