function [ filt info ] = globalmax_filt(ds, orig, params)
% brightest point - preserves brightest point's value in filt, returns
% point coord in info.center array of Nx3 (ZYX coord for each syn)
if(strcmp(params.type, 'globalmax'))
    info.center = zeros(ds.ntrain,3); % Z Y X coords of brightest point
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain
        slice = orig(s,:,:,:);
        [val ind] = max(slice(:));
        [z y x] = ind2sub(ds.sgdim, ind);
        info.center(s,:) = [z y x];
        filt(s,z,y,x) = val;
    end
    return;
end


