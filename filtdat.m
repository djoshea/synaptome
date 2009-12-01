function [filt info] = filtdat(ds, dat, params)
% returns a filtered version of data
% dat is either a channel name (in trdlist) or a ntrain x sgdim... array
% params is a struct with field type and other type-specific settings
% filt is an array ntrain x sgdim1 x ...

info = [];
filt = [];

if(~exist('params', 'var') || ~isfield(params, 'type'))
    error('Filter type not specified in params.type');
end

if(ischar(dat))
    % get data channel by name
    orig = getchannel(ds, dat);
else
    % assume dat is the data
    orig = dat;
end

% ordering is due to mask(z, y, x) indexing
[Y Z X] = meshgrid(1:ds.sgdim(1), 1:ds.sgdim(2), 1:ds.sgdim(3));
distmask = @(center) sqrt(((Z-center(1))*ds.sgaspect(1)).^2 + ...
                          ((Y-center(2))*ds.sgaspect(2)).^2 + ...
                          ((X-center(3))*ds.sgaspect(3)).^2 ); 

% simple binary masking
if(strcmp(params.type, 'mask'))
    if(isfield(params, 'mask'))
        mask = params.mask;
        if(ndims(mask) == 2)
            mask = repmat(shiftdim(mask,-2), [ds.ntrain ds.sgdim(1) 1 1]);
        end
        if(ndims(mask) == 3)
            mask = repmat(shiftdim(mask,-1), [ds.ntrain 1 1 1]);
        end
    else
        error('No params.mask specified for filter mask.');
    end
    
    filt = mask .* orig;
    return;
end
                      
% Mask ball filtering
% center and radius can be a 1x3 array and scalar for same center in all
% synapses, or they can be an Nx3 and
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

% Tophat filtering
if(strcmp(params.type, 'tophat'))
    if(isfield(params, 'se'))
        se = params.se;
    elseif(isfield(params, 'pxradius3'))
        mask = distmask(center) <= params.pxradius3;
        se = strel(mask);
    elseif(isfield(params, 'pxradius'))
        se = strel('disk', params.pxradius, 0);
    else
        pxradius = 2;
        se = strel('disk', pxradius, 0);
    end
    
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain    
        % filter each synapse independently
        % shift dim goes from ZYX indexing to YXZ indexing and back
        sh = shiftdim(squeeze(orig(s,:,:,:)),1);
        filt(s,:,:,:) = shiftdim(imtophat(sh, se), 2);
    end
    return;
end

% Maxima
if(strcmp(params.type, 'regmax'))    
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain    
        % filter each synapse independently
        % shift dim goes from ZYX indexing to YXZ indexing and back
        sh = shiftdim(squeeze(orig(s,:,:,:)),1);
        filt(s,:,:,:) = shiftdim(imregionalmax(sh), 2);
    end
    return;
end

% h-Maxima
if(strcmp(params.type, 'hmax'))    
    if(isfield(params, 'h'))
        h = params.h;
    else
        h = 0.2;
    end
    
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain    
        % filter each synapse independently
        % shift dim goes from ZYX indexing to YXZ indexing and back
        sh = shiftdim(squeeze(orig(s,:,:,:)),1);
        filt(s,:,:,:) = shiftdim(imhmax(sh,h), 2);
    end
    return;
end

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

% Maxima within neighborhood 2D only
if(strcmp(params.type, 'maxnhood'))
    if(isfield(params, 'se'))
        se = params.se;
    elseif(isfield(params, 'pxradius'))
        se = strel('disk', params.pxradius, 0);
    elseif(isfield(params, 'radius'))
        se = strel('disk', ceil(params.radius / ds.sgaspect(1)), 0);
    end
    
    nhood = getnhood(se);
    order = nnz(nhood);
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain
        for z = 1:ds.sgdim(1)
            slice = squeeze(orig(s,z,:,:));
            filt(s,z,:,:) = slice .* (slice == ordfilt2(slice, order, nhood));
        end
    end
    return;
end


% Maxima within 3D Window
if(strcmp(params.type, 'maxwind3'))
    if(isfield(params, 'pxradius'))
        pxradius = params.pxradius;
    else
        pxradius = 1;
    end
    if(isfield(params, 'thresh'))
        thresh = params.thresh;
    else
        thresh = 0;
    end
    
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain
        slice = squeeze(orig(s,:,:,:));
        filt(s,:,:,:) = (slice == ordfilt3(slice, 'max', pxradius)).* slice >= thresh ;
    end
    return;
end

% Voronoi mask: return mask that selects pixels closer to me 
% (pivot closest to centermost pixel) than any other pixel in
% the input channel (should be binary)
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

% Watershed mask: return mask that selects pixels within the catchment basin 
% belonging to the marker closest to the params.center (or center of volume).
% params.markers is ds.ntrain x ds.sgdim binary array with markers at the 
% maxima to watershed from. 
% input channel is the grayscale channel to watershed upon.
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



error(sprintf('Filter "%s" not implemented', params.type));