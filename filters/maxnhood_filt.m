function [ filt info ] = maxnhood_filt( ds, orig, params )

info = [];
filt = [];

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
end

