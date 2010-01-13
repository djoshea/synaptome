syn = getchannel(ds, 'Synapsin', 119);
synmax = getchannel(ds, 'Synapsin_mw3', 119);

syninv = -syn + max(syn(:));
synimposed = imimposemin(syninv, synmax);

L = watershed(synimposed);

Lflat = reshape(L,size(L,1),[]);
RGB = ind2rgb(Lflat, [0 0 0; jet(length(unique(L))-1)]);

C = zeros([size(L) 3]);
for c = 1:3
    C(:,:,:,c) = reshape(RGB(:,:,c), size(L));
end

viewsgrow(C);