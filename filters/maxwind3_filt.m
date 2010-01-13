function [ filt info] = maxwind3_filt( ds, orig, params )

filt = [];
info = [];

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


end

