function [ filt info ] = maxseed_filt( ds, orig, params )
% Maxima that exceed threshold (auto determined by graythresh or
% params.thresh) that are far enough away to not be fused by opening with
% a se (params.se or cube of size 3x3x3, i.e. unique within 5x5x5
% neighborhood)

filt = [];
info = [];

if(strcmp(params.type, 'maxseed'))
    if(~isfield(params,'thresh'))
        threshauto = 1;
    else
        threshauto = 0;
    end
    if(~isfield(params, 'se'))
        se = strel('arbitrary', ones(3,3,3));
    else
        se = params.se;
    end
    
    filt = zeros([ds.ntrain ds.sgdim]);
    for s = 1:ds.ntrain    
       syn = squeeze(orig(s,:,:,:));
       if(threshauto)
           thresh = graythresh(syn);
       else
           thresh = params.thresh;
       end
       seeds = imregionalmax(syn) .* (syn > thresh);
       dil = imdilate(seeds, se);
       prop = regionprops(bwlabeln(dil),seeds,'WeightedCentroid');
       centroid = round(cat(1,prop.WeightedCentroid));
       inds = sub2ind(ds.sgdim, centroid(:,2), centroid(:,1), centroid(:,3));
       
       centroids = zeros(ds.sgdim);
       centroids(inds) = 1;
       filt(s,:,:,:) = centroids;
    end
    return;
end

