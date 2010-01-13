function [ filt info ] = mask_filt(ds, orig, params)
% Simple binary masking

filt = [];
info = [];

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

end

