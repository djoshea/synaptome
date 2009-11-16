function sg = loadsg(ds, i)
% loadsg(ds, i) - Load synaptogram i from ds.trainPivotss

if(i < 1 || i > length(ds.trainPivots.x))
    error('Invalid synapse index');
end

sg.i = i;
sg.x = ds.trainPivots.x(i);
sg.y = ds.trainPivots.y(i);
sg.z = ds.trainPivots.z(i);

ds.sgdim = [11 11 11]; % x, y ,z
ds.sgdir = 'U:\Brad\gordon 20090630 32 bit\tmp\';
% ds.sgdir = '/Volumes/experiments/Brad/gordon 20090630 32 bit/tmp/';

sg.im = zeros(ds.nch, ds.sgdim(3), ds.sgdim(2), ds.sgdim(1)); % ch, z, y, x

width = max(4, max(ceil(log10(double([sg.x, sg.y, sg.z]) + 1))));
fmatcoords = sprintf('%%0%dd-%%0%dd-%%0%dd', width, width, width);
sg.coordstr = sprintf(fmatcoords, sg.z, sg.y, sg.x);

for c = 1:ds.nch
    fname = sprintf('%s%s.tif %s.tif', ds.sgdir, ds.chlist{c}, sg.coordstr);
    
    for z = 1:ds.sgdim(3)
        sg.im(c,z,:,:) = imread(fname, z);
    end
    
    
end