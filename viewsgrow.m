function rowdat = viewsgrow(dat, bordercol, vspacing, hspacing)
% dat given in (z,y,x) coords, or specify dat.R .G .B separately

if(~isstruct(dat)) % just data, make into matrix
    d.R = squeeze(dat);
    d.G = squeeze(dat);
    d.B = squeeze(dat);
else
    fn = fieldnames(dat);
    z = zeros(size( squeeze(getfield(dat,char(fn(1)))) ));
    if(isfield(dat, 'R'))
        d.R = squeeze(dat.R);
    else
        d.R = z;
    end
    if(isfield(dat, 'G'))
        d.G = squeeze(dat.G);
    else
        d.G = z;
    end
    if(isfield(dat, 'B'))
        d.B = squeeze(dat.B);
    else
        d.B = z;
    end
end

ntiles = size(d.R,1); % Z
szy = size(d.R,2);
szx = size(d.R,3);

if(~exist('bordercol', 'var'))
    bordercol = [0 0 0];
end
if(~exist('vspacing', 'var'))
    vspacing = 1;
end
if(~exist('hspacing', 'var'))
    hspacing = 2*vspacing;
end

% fill rowdat (Y, X, C) with bordercol
rowdat = repmat(reshape(bordercol, [1 1 3]), [2*vspacing + szy, ntiles*(szx+hspacing) + hspacing, 1]);

% should be 1 if channels are normalized;
norm = max([d.R(:); d.G(:); d.B(:)]);
for z = 1:ntiles
    yind = vspacing+1 : vspacing+szy;
    xind = szx*(z-1)+hspacing*z+1 : szx*z+hspacing*z;
%     rowdat(yind, xind, 1) = squeeze(d.R(z,:,:)) / max(d.R(:));
%     rowdat(yind, xind, 2) = squeeze(d.G(z,:,:)) / max(d.G(:));
%     rowdat(yind, xind, 3) = squeeze(d.B(z,:,:)) / max(d.B(:));
    rowdat(yind, xind, 1) = squeeze(d.R(z,:,:)) / norm;
    rowdat(yind, xind, 2) = squeeze(d.G(z,:,:)) / norm;
    rowdat(yind, xind, 3) = squeeze(d.B(z,:,:)) / norm;
end

rowdat = min(max(rowdat, 0), 1);

if(nargout == 0)
    figure;
    sh = rowdat - min(rowdat(:));
    sh = sh / max(sh(:));
    imagesc(sh);
end
    
