function dat = getchannel(ds, dname, syn)
% Retrieves a specific data channel by name from ds.trd array
% if syn is specified, gets a certain synapse's data
% chdat will be an array: ntrain x sgdim(1) x sgdim(2) x sgdim(3) [ x 3 for color]

if(~strcmp(dname, ''))
    chidx = find(strcmp(ds.trdlist, dname));
    if(numel(chidx) == 0)
        % couldn't find in trdlist, check im channels?
        dat = getimchannel(ds, dname);
        if(isnan(dat))
            % try color channels
            chidx = find(strcmp(ds.colorchname, dname));
            if(numel(chidx) == 0)
                % still nothing
                error('Could not find channel %s', dname);
            end
            dat = squeeze(ds.colorch(chidx(1),:,:,:,:,:)); % grab the color channel
        end
    else
        dat = squeeze(ds.trd(chidx,:,:,:,:)); % grab the non-color channel
    end
else
    % fill with zeros if '' is specified for the name
    dat = zeros([ds.ntrain ds.sgdim]);
end
 
% grab a specific synapse if requested    
if(exist('syn', 'var'))
    dat = squeeze(dat(syn,:,:,:,:));
end


