function ds = addcolorchannel(ds, dat, dname)
% Adds a channel to the ds.colorch array
% dat should be an array: ntrain x 11 x 11 x 11 x 3

if(~isfield(ds, 'colorch'))
    ds.colorch = zeros([0 ds.ntrain ds.sgdim 3]);
    ds.colorchname = {};
end

ds.colorch(end+1,:,:,:,:,:) = shiftdim(dat, -1);
ds.colorchname{end+1} = char(dname);