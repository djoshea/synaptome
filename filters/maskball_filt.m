function [ filt info ] = maskball_filt(ds, orig, params)
%MASKBALL_FILT Mask ball filtering
% center and radius can be a 1x3 array and scalar for same center in all
% synapses, or they can be an Nx3 to specify different centers and radii

filt = [];
info = [];

[Y Z X] = meshgrid(1:ds.sgdim(1), 1:ds.sgdim(2), 1:ds.sgdim(3));
distmask = @(center) sqrt(((Z-center(1))*ds.sgaspect(1)).^2 + ...
                          ((Y-center(2))*ds.sgaspect(2)).^2 + ...
                          ((X-center(3))*ds.sgaspect(3)).^2 );                      

if(strcmp(params.type, 'maskball'))
    if(isfield(params, 'center'))
        center = params.center;
    else
        % default is center of volume
        center = ds.sgdim/2 + 1/2;
    end
    
    if(isfield(params, 'radius'))
        radius = params.radius; %should be in same units as ds.sgaspect
    else
        % assume they meant largest sphere that fits in synaptogram
        radius = min((ds.sgaspect.*ds.sgdim)/2 - 1/2);
    end
    
    if(size(center,1) ~= ds.ntrain)
        center = repmat(center, [ds.ntrain 1]);
    end
    if(numel(radius) ~= ds.ntrain)
        radius = repmat(radius, [ds.ntrain 1]);
    end

    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain
        dmask = distmask(center(s,:));
        mask = distmask(center(s,:)) <= radius(s);
        filt(s,:,:,:) = squeeze(orig(s,:,:,:)) .* mask;
    end
    return;
end


