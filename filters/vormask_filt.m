function [ filt info ] = vormask_filt( ds, orig, params)
% Voronoi mask: return mask that selects pixels closer to me 
% (pivot closest to centermost pixel) than any other pixel in
% the input channel (should be binary)

filt = [];
info = [];

[Y Z X] = meshgrid(1:ds.sgdim(1), 1:ds.sgdim(2), 1:ds.sgdim(3));
distmask = @(center) sqrt(((Z-center(1))*ds.sgaspect(1)).^2 + ...
                          ((Y-center(2))*ds.sgaspect(2)).^2 + ...
                          ((X-center(3))*ds.sgaspect(3)).^2 ); 

if(strcmp(params.type, 'vormask'))
    if(isfield(params, 'biasfactor'))
        biasfactor = params.biasfactor;
    else
        biasfactor = 1;
    end
    
    % find pivot closest to center
    center = ds.sgdim/2 + 1/2;
    dmask = distmask(center);
 
    % for storing coords of pivot closest to centroid
    info.center = zeros(ds.ntrain, length(ds.sgdim));
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain
        slice = squeeze(orig(s,:,:,:)) > 0;
        
        mm = (2*max(dmask(:)) - dmask) .* slice;
        [val ind] = max(mm(:));
        [z y x] = ind2sub(ds.sgdim, ind);
        info.center(s,:) = [z y x];
        info.distcenter(s,:) = dmask(z,y,x);
        
        % pivots sans central pivot
        others = slice;
        others(ind) = 0;
        othersdist = bwdist(others);
        
        central = zeros(size(slice));
        central(ind) = 1;
        centraldist = bwdist(central);
        
         filt(s,:,:,:) = othersdist > biasfactor*centraldist;
    end

    return;
end
