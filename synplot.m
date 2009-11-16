function synplot(ds, features)

features = {'IB_VGlut1', 'IB_PSD95'};

[ft idx] = getfeature(ds, features);
    
clf;
set(gcf, 'defaulttextinterpreter', 'none');

for ci = 1:ds.nlabels
    pts = ds.trainlabel == ci;
    for i = 1:nnz(pts)
        idx = pts(i);
        h = line(ft(idx, 1), ft(idx, 2));
        set(h, 'Marker','.', 'LineStyle', 'none', ...
         'Color', ds.labelcolors(ci,:), ...
         'MarkerFaceColor', ds.labelcolors(ci,:), ...
         'MarkerEdgeColor', ds.labelcolors(ci,:), ...
         'MarkerSize', 8, ...
         'ButtonDownFcn', @callback);
    hold on
end

xlabel(features{1});
ylabel(features{2});
title(sprintf('%s vs. %s', features{1}, features{2}));

function callback(gcbo, eventdata, handles)

srcdisp('hello!')