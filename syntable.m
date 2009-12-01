function syntable(ds, i)

% data table figure with synapse features
ti = sprintf('Data for Synapse %d [ %s ]', i, ds.sg(i).str);

figure(198); clf;
% set(gcf, 'Position', [1661 476 350 400]);
set(gcf, 'Name', ti);
set(gcf, 'NumberTitle', 'off');

t = uitable;
set(t, 'Units', 'normalized')
set(t, 'Position', [0.01 0.01 0.98 0.98]);

set(t,'RowStriping', 'on');
set(t,'BackgroundColor',[1 1 1; 0.9 0.9 1]);

set(t,'ColumnName', {'Feature', 'Value', 'Class Mean', 'Overall Mean'});
ftdat = mat2cell(ds.ft(i,:)', ones(size(ds.ft,2), 1));
allmean = mean(ds.ft, 1)';
classmean = mean(ds.ft(ds.trainlabel == ds.trainlabel(i),:), 1)'; 

dat = [ds.ftname' ftdat ...
    mat2cell(classmean, ones(size(classmean))) ...
    mat2cell(allmean, ones(size(allmean)))];
set(t,'Data',dat);