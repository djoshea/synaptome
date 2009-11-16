function pivots = loadpivotlist(fnamepivot)

% fnamepivot = '/Volumes/experiments/Brad/gordon 20090630 32 bit/Synapsin.tif.pivots.xls';
fpivot = fopen(fnamepivot);
data = textscan(fpivot,'%d %*d %f %d %d %d', 'headerlines',1);
fclose(fpivot);

pivots.idx = data{1};
pivots.brightness = data{2};
pivots.x = data{3};
pivots.y = data{4};
pivots.z = data{5};

