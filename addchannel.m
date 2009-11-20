function ds = addchannel(ds, dat, dname)
% Adds a channel to the ds.trdat array
% dat should be an array: ntrain x 11 x 11 x 11

ds.trd(end+1,:,:,:,:) = shiftdim(dat, -1);
ds.trdlist{end+1} = char(dname);