%% Create data channels from imaging data
% this involves copying processed data from ds.sg.im to ds.dat
% use addchannel(ds, chdat, chname) to add to ds.trd

% ds.ch is a ndatachannels x ntrain x sgdim 5-D array of processed data
ds.trd = zeros([0 ds.ntrain ds.sgdim]);
ds.trdlist = cell(0);

% for now, just copy raw image data as channels
pfilt = [];
for i = 1:length(ds.chlist)
    chname = ds.chlist{i};
    ds = addchannel(ds, getchannel(ds, chname), chname);
end

% top hat filter synapsin
ftophat = struct('type', 'tophat', 'se', strel('disk',2));
ds = addchannel(ds, filtdat(ds, 'Synapsin', ftophat), 'Synapsin_th');

% then find maxima in 5x5x5 neighborhoods
fmaxwind3 = struct('type', 'maxwind3', 'pxradius', 5, 'thresh', 0.05);
ds = addchannel(ds, filtdat(ds, 'Synapsin', fmaxwind3), 'Synapsin_mw3');

% construct a voronoi region around the central-most synapsin puncta
% store central synapsin puncta coords in info.center
% this mask avoids other synapsin puncta
fwatershed = struct('type', 'watershedmask', 'markers', getchannel(ds, 'Synapsin_mw3'));
[watermask info] = filtdat(ds, 'Synapsin', fwatershed);
synapsin_cent = info.centermarker;
synapsin_dist = info.distcentermarker;
fwatermask = struct('type', 'mask', 'mask', watermask);

% ball restriction mask around central synapsin puncta, further restricts
% voronoi mask that avoids other synapsin puncta
frestrictpre = struct('type', 'maskball', 'radius', 250, 'center', synapsin_cent);
frestrictpost = struct('type', 'maskball', 'radius', 400, 'center', synapsin_cent);

% Final masks for pre and post, ball restricted and voronoi restricted
maskpre = filtdat(ds, watermask, frestrictpre);
fmaskpre = struct('type', 'mask', 'mask', maskpre);
maskpost = filtdat(ds, watermask, frestrictpost);
fmaskpost = struct('type', 'mask', 'mask', maskpost);

% to show final masks over synapsin channel
ds = addchannel(ds, watermask, 'Synapsin_watermask');
ds = addchannel(ds, filtdat(ds, 'Synapsin', fmaskpre), 'Synapsin_maskpre');
ds = addchannel(ds, filtdat(ds, 'Synapsin', fmaskpost), 'Synapsin_maskpost');

% global max for brightest puncta find
fmax = struct('type', 'globalmax');
% small ball around global max to pass only active brightness zone
fab = struct('type', 'maskball', 'radius', 200);

% active brightness filtering: 
% search for brightest point (fmask) within masked region (fmaskpre or fmaskpost)
% pass only a small ball (fab) around the brightest point
% later this active brightness zone will be integrated to form the _iab features
ablistpre = {'Bassoon','VGlut1','VGlut2','VGat', 'GAD'};
ablistpost = { 'PSD95', 'Gephyrin' };
ablist = [ablistpre ablistpost];

centdists = zeros(ds.ntrain, length(ablist));
distfn = @(zyx,center) sqrt(((zyx(:,1)-center(:,1))*ds.sgaspect(1)).^2 + ...
                          ((zyx(:,2)-center(:,2))*ds.sgaspect(2)).^2 + ...
                          ((zyx(:,3)-center(:,3))*ds.sgaspect(3)).^2 );
syncenter = ds.sgdim/2 + 1/2; 
    
fprintf('Active Brightness: ');
for c = 1:length(ablist)
    cname = ablist{c};
    fprintf('%s ', cname);
    cname_th = sprintf('%s_th',cname);
    cname_ab = sprintf('%s_ab',cname);
    
%     ds = addchannel(ds, filtdat(ds, cname, ftophat), cname_th);
    
    % find all maxima using neighborhood maxima
    maxima = filtdat(ds, cname, fmaxwind3);
    
    % find watershed basin mask for maximum closest to synapsin centroid
    fwater = struct('type', 'watershedmask', 'center', synapsin_cent, ...
        'markers', maxima);
    [filt info] = filtdat(ds, cname, fwater);
    
    % distance to maximum closest to synapsin centroid, used as feature
    centdists(:, c) = distfn(info.centermarker, synapsin_cent);
    
    if(nnz(strcmp(ablistpre, cname))) % use pre or post synaptic search restriction mask?
         restricted = filtdat(ds, filt, fmaskpre);
    else
         restricted = filtdat(ds, filt, fmaskpost);
    end
    
    fmask = struct('type', 'mask', 'mask', restricted);
    masked = filtdat(ds, cname, fmask);
    ds = addchannel(ds, masked, cname_ab);
end
fprintf('\n');

%% Specify synaptogram visualization structure
% this involves creating a cell array with one element for each row
% each row consists of 1 or 3 channel names that index into ds.ch
% 1 name implies grayscale, 3 names are used as RGB channels
ds.vis = {};
ds.visname = {};

vis = { 'Synapsin', 'Synapsin_mw3', 'Synapsin_watermask', 'Synapsin_maskpre', 'Synapsin_maskpost', ...
       'Bassoon', ...
    'VGlut1', {'Synapsin_maskpre','VGlut1_ab'}, ...
    'VGlut2', {'Synapsin_maskpre', 'VGlut2_ab'}, ...
    'PSD95', {'Synapsin_maskpost', '', 'PSD95_ab'}, ...
    'VGat', {'Synapsin_maskpre','VGat_ab'}, ...
    'GAD', {'Synapsin_maskpre', 'GAD_ab'}, ...
    'Gephyrin', {'Synapsin_maskpost', '', 'Gephyrin_ab'}, ...
    };

for i = 1:length(vis)
    ds = addvisrow(ds, vis{i});
end

ds.ft = [];
ds.ftname = {};

% centroid distance features
ds = addfeature(ds, synapsin_dist_center, 'Synapsin_dist');
for c = 1:length(ablist)
    ftname = sprintf('%s_dist', ablist{c});
    ds = addfeature(ds, centdists(:,c), ftname);
end

% integrated brightness features
% for c = 1:ds.nimch
%    cname = ds.chlist{c};
%    name = sprintf('%s_ib', cname);
%    dat = zeros(ds.ntrain,1);
%    for i = 1:ds.ntrain
%        chdat = min(ds.sg(i).im(6,5:7,5:7,5:7), ds.sg(i).im(c,5:7,5:7,5:7));
%        dat(i) = sum(chdat(:));
%    end
%    
%    ds = addfeature(ds, dat, name);
% end

% active brightness features
for c = 1:length(ablist)
   cname = ablist{c};
   name = sprintf('%s_iab', cname);
   dat = getchannel(ds, sprintf('%s_ab', cname));
   ftdat = zeros(ds.ntrain,1);
   for s = 1:ds.ntrain
       dslice = dat(s,:,:,:);
       ftdat(s) = sum(dslice(:));
   end
   ds = addfeature(ds, ftdat, name);
end

glutdat = getfeature(ds,{'PSD95_iab','VGlut1_iab','VGlut2_iab'});
glutdat = glutdat - repmat(min(glutdat,[],1),[ds.ntrain 1]);
glutdat = glutdat ./ repmat(max(glutdat,[],1),[ds.ntrain 1]);
gabadat = getfeature(ds,{'Gephyrin_iab','VGat_iab', 'GAD_iab'});
gabadat = gabadat - repmat(min(gabadat,[],1),[ds.ntrain 1]);
gabadat = gabadat ./ repmat(max(gabadat,[],1),[ds.ntrain 1]);

ds = addfeature(ds, sqrt(sum(glutdat(:,2:3).^2,2)), 'Glut_L2pre');
ds = addfeature(ds, sqrt(sum(gabadat(:,2:3).^2,2)), 'GABA_L2pre');


ds = addfeature(ds, min([glutdat(:,1) sqrt(sum(glutdat(:,2:3).^2,2))], [], 2), 'Glut_L2prepost');
ds = addfeature(ds, min([gabadat(:,1) sqrt(sum(gabadat(:,2:3).^2,2))], [], 2), 'GABA_L2prepost');


% which training examples to actually plot and test on
ds.trainactive = ds.trainlabelconf > 0.4;
% ds.trainactive = ds.trainactive & ds.trainlabel ~= 1 & ds.trainlabel ~= 4;

%% Test CV error
% testft;

%% Plot a 2D separability scatter plot
figure(1), clf;
feat = {'GABA_L2pre', 'Gephyrin_iab'};
synplot(ds, feat);

figure(2), clf;
feat = {'Glut_L2pre', 'PSD95_iab'};
synplot(ds, feat);
