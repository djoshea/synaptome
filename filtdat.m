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
        se = strel('disk', params.pxradius);
    else
        pxradius = 2;
        se = strel('disk', pxradius);
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

error(sprintf('Filter "%s" not implemented', params.type));