function chdat = getimchannel(ds, chname)
% Retrieves a specific imaging channel by name from ds.sg(:).im structure
% chdat will be an array: ntrain x sgdim(1) x sgdim(2) x sgdim(3) ...

chdat = zeros([ds.ntrain ds.sgdim]);
chidx = find(strcmp(ds.chlist, chname));
if(numel(chidx) == 0)
    chdat = NaN;
    return;
end

for i = 1:ds.ntrain
    chdat(i,:,:,:) = ds.sg(i).im(chidx, :,:,:);
end


