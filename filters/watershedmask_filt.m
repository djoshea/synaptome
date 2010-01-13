function [ filt info ] = watershedmask_filt(ds, orig, params )
% Watershed mask: return mask that selects pixels within the catchment basin 
% belonging to the marker closest to the params.center (or center of volume).
% params.markers is ds.ntrain x ds.sgdim binary array with markers at the 
% maxima to watershed from. 
% input channel is the grayscale channel to watershed upon.

filt = [];
info = [];

if(strcmp(params.type, 'watershedmask'))
    % find pivot closest to center
    if(~isfield(params, 'center'))
        center = repmat(ds.sgdim/2 + 1/2, [ds.ntrain 1]); % choose center of volume
    else
        if(ndims(params.center) == 1) % clone by ntrain
            center = repmat(params.center, [ds.ntrain 1]);
        else
            center = params.center;
        end
    end
    
    if(~isfield(params, 'markers'))
        error('Must specify markers in params.markers array');
    end
 
    % for storing coords of pivot closest to centroid
    info.centermarker = zeros(ds.ntrain, length(ds.sgdim));
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain
        slice = squeeze(orig(s,:,:,:));
        markerslice = squeeze(params.markers(s,:,:,:) > 0);
        
        % compute marker closest to center and distance to center
        dmask = distmask(center(s,:));
        mm = (2*max(dmask(:)) - dmask) .* markerslice;
        [val ind] = max(mm(:));
        [z y x] = ind2sub(ds.sgdim, ind);
        info.centermarker(s,:) = [z y x];
        info.distcentermarker(s,:) = dmask(z,y,x);
        
        sliceinv = -slice + max(slice(:));
        sliceimposed = imimposemin(sliceinv, markerslice);
        L = watershed(sliceimposed);
        filt(s,:,:,:) = L==L(z,y,x);
    end

    return;
end

