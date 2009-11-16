% data table figure with synapse features
i = 1;
ti = sprintf('Data for Synapse %d [ %s ]', i, sg.coordstr);

figure(198); clf;
set(gcf, 'Position', [1661 476 450 700]);
set(gcf, 'Name', ti);
set(gcf, 'NumberTitle', 'off');

t = uitable;
set(t, 'Units', 'normalized')
set(t, 'Position', [0.01 0.01 0.98 0.98]);

ftdat = mat2cell(ds.ft(i,:)', ones(size(ds.ft,2), 1));
dat = {ds.ftname, ftdat};
set(t,'Data',dat);