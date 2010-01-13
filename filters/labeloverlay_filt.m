function [ filt info ] = labeloverlay_filt( ds, orig, params)
% Color a grayscale channel in regions marked by params.label matrix
% specify NumLabels (excluding 0) x 3 RGB colormap in params.cmap, default is hsv(NumLabels)
% with 1 1 1 (original grayscale) for label=0
% filt is ds.ntrain x ds.sgdim x 3 (5-D)

info = [];
filt = [];

if(strcmp(params.type, 'labeloverlay'))
    if(~isfield(params, 'label'))
        error('Must specify label matrix in params.label');
    else
        L = params.label;
    end
    
    info.cmap = cell(0);
    filt = repmat(orig,[ones(1,ndims(orig)) 3]);
    for s=1:ds.ntrain
        NL = max(L(s,:)); % number 
        if(~isfield(params, 'cmap'))
            cmap = jet(NL);
            cmap = cmap(randperm(NL),:);
        elseif(strcmp(params.cmap, 'function_handle'))
            cmap = params.cmap(NL);
            cmap = cmap(randperm(NL),:);
        else
            cmap = params.cmap{s};
            if(size(cmap,1) < NL)
                error('Not enough colors in colormap');
            end
        end
        
        info.cmap{s} = cmap;
        
        for z = 1:ds.sgdim(1) % loop over 2D images to be able to use built-ins
            slice = squeeze(orig(s,z,:,:)); % 2D xy image
            Lslice = squeeze(L(s,z,:,:));
            
            % scale the brightness value by the original image grayscale
            % keep the hue and saturation the same as in cmap
            cslicehsv = rgb2hsv(label2rgb(Lslice,cmap,[1 1 1]));
            cslicehsv(:,:,3) = cslicehsv(:,:,3) .* slice;
            filt(s,z,:,:,:) = hsv2rgb(cslicehsv);
        end
    end
    return;
end

% extended maxima (regional maxima of hmaxima transform)
% then dilates then erodes with 2x2x2 cube to link adjacent maxima
if(strcmp(params.type, 'hmax'))
    if(isfield(params, 'h'))
        h = params.h;
    else
        error('Must specify params.h');
    end
    
    se = strel('arbitrary', ones(2,2,2));
    filt = zeros(size(orig));
    for s= 1:ds.ntrain
        allmax = imextendedmax(squeeze(orig(s,:,:,:)), h);
        filt(s,:,:,:) = imerode(imdilate(allmax,se),se);
    end
    return;
end
