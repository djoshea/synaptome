function [ filt info ] = hmax_filt(ds, orig, params)
% h-Maxima: see 'help imhmax' with params.h as h value

filt = [];
info = [];

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

