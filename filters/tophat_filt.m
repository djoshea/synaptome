function [ filt info ] = tophat_filt( ds, orig, params )
% Tophat filtering

filt = [];
info = [];

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


