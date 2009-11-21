function ds = addchannel(ds, dat, dname, normalize)
% Adds a channel to the ds.trdat array
% dat should be an array: ntrain x 11 x 11 x 11
% normalize (default on) will shift linearly into [0,1] range

if(~exist('normalize', 'var'))
    normalize = 1;
end
if(normalize)
    dat = dat - min(dat(:));
    dat = dat / max(dat(:));
end

ds.trd(end+1,:,:,:,:) = shiftdim(dat, -1);
ds.trdlist{end+1} = char(dname);