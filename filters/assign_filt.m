function [ filt info ] = assign_filt(ds, orig, params )
% Assign: Takes original channel and params.sourceRegions label matrix, computes the
% weighted centroid of each connected component (puncta), samples the 
% params.targetRegions label matrix to compute the region to which each centroid
% belongs, and relabels the original channel with the label ids of the
% params.regions matrix to which each blob's centroid belongs

filt = [];
info = [];

if(strcmp(params.type, 'assign'))
    if(~isfield(params,'sourceRegions'))
        error('Must specify params.sourceRegions label matrix');
    end
    if(~isfield(params,'targetRegions'))
        error('Must specify params.targetRegions label matrix');
    end
    
    filt = zeros([ds.ntrain ds.sgdim]);
    for s=1:ds.ntrain
        sourceSlice = squeeze(params.sourceRegions(s,:,:,:));
        targetSlice = squeeze(params.targetRegions(s,:,:,:));
        props = regionprops(sourceSlice, squeeze(orig(s,:,:,:)), 'WeightedCentroid');
        centroids = round(cat(1, props.WeightedCentroid));
        % grab the target label id for each source label id
        Lcentroids = targetSlice(sub2ind(ds.sgdim, centroids(:,2), centroids(:,1), centroids(:,3)));
        filtslice = zeros(ds.sgdim);
        for val = 1:size(centroids, 1) % replace each instance of val in source slice with the target centroid value
           filtslice(sourceSlice == val) = Lcentroids(val);
        end
        filt(s,:,:,:) = filtslice;
    end
    return
end