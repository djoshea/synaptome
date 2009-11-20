function dat = getchannel(ds, dname, syn)
% Retrieves a specific data channel by name from ds.ch array
% if syn is specified, gets a certain synapse's data
% chdat will be an array: ntrain x sgdim(1) x sgdim(2) x sgdim(3) ...

if(~strcmp(dname, ''))
    chidx = find(strcmp(ds.trdlist, dname));
    if(numel(chidx) == 0)
        error(sprintf('Could not find channel "%s"', dname));
    end
    dat = squeeze(ds.trd(chidx,:,:,:,:));
else
    % fill with zeros if '' is specified for the name
    dat = zeros([ds.ntrain ds.sgdim]);
end
 
% grab a specific synapse if requested    
if(exist('syn', 'var'))
    dat = dat(syn,:,:,:);
end


