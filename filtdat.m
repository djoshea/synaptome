function filt = filtdat(ds, dname, params)
% returns a filtered version of channel dname
% params is a struct with field type and other type-specific settings
% filt is an array ntrain x sgdim1 x ...

orig = getchannel(ds, dname);

if(strcmp(params.type, 'maskball'))
    center = params.center;
    rad = params.rad;
    
    