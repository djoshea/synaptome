function [ filt info ] = watershedlabel_filt(ds, orig, params)
% Watershed label: return label matrix marking each catchment basin 
% params.markers is ds.ntrain x ds.sgdim binary array with markers at the 
% maxima to watershed from. 
% input channel is the grayscale channel to watershed upon.

filt = [];
info = [];

[Y Z X] = meshgrid(1:ds.sgdim(1), 1:ds.sgdim(2), 1:ds.sgdim(3));
distmask = @(center) sqrt(((Z-center(1))*ds.sgaspect(1)).^2 + ...
                          ((Y-center(2))*ds.sgaspect(2)).^2 + ...
                          ((X-center(3))*ds.sgaspect(3)).^2 );     

if(strcmp(params.type, 'watershedlabel'))
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
        % supersample to thin boundary between catchment basins (to no
        % width)
%         N = 3;
%         indZ = 1:ds.sgdim(1);
%         indZ = reshape(indZ(ones(N,1),:),1,[]);
%         indY = 1:ds.sgdim(2);
%         indY = reshape(indY(ones(N,1),:),1,[]);
%         indX = 1:ds.sgdim(3);
%         indX = reshape(indX(ones(N,1),:),1,[]);
%         [Z Y X] = meshgrid(indZ,indY,indX);
%         slicebig = sliceimposed(sub2ind(ds.sgdim,Z,Y,X));
%         
%         wshed = watershed(slicebig);
%         L = wshed(ceil(N/2):N:end, ceil(N/2):N:end, ceil(N/2):N:end);
        L = watershed(sliceimposed);
        [D Ind] = bwdist(L>0);
        centerlabel = L(z,y,x);
        
        % assign label = 1 to region closest to synapsin, swap with other
        % label
        L = L(Ind);
        Lcopy = L;
        L(Lcopy == centerlabel) = 1; 
        L(Lcopy == 1) = centerlabel;
        info.centerlabel = 1;
        
        filt(s,:,:,:) = L;
    end

    return;
end