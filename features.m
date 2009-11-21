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



% active brightness filtering
ftophat = struct('type', 'tophat', 'se', strel('disk',2));
frestrict = struct('type', 'maskball', 'radius', 300);
fmax = struct('type', 'globalmax');
fab = struct('type', 'maskball', 'radius', 200);

ds = addchannel(ds, filtdat(ds, 'Synapsin', frestrict), 'Synapsin_rs');

ablist = {'VGlut1','VGlut2', 'PSD95', 'VGat', 'GAD', 'Gephyrin'};

for c = 1:length(ablist)
    cname = ablist{c};
    cname_th = sprintf('%s_th',cname);
    cname_ab = sprintf('%s_ab',cname);
    
    ds = addchannel(ds, filtdat(ds, cname, ftophat), cname_th);
    restricted = filtdat(ds, cname_th, frestrict);
    [filt info] = filtdat(ds, restricted, fmax);
    fab.center = info.center;
    ds = addchannel(ds, filtdat(ds, cname_th, fab), cname_ab);
end

%% Specify synaptogram visualization structure
% this involves creating a cell array with one element for each row
% each row consists of 1 or 3 channel names that index into ds.ch
% 1 name implies grayscale, 3 names are used as RGB channels

ds.vis = {};
ds.visname = {};

vis = { 'Synapsin', 'Bassoon', ...
    'VGlut1', {'Synapsin_rs','VGlut1_ab'}, ...
    'VGlut2', {'Synapsin_rs', 'VGlut2_ab'}, ...
    'PSD95', {'Synapsin_rs', 'PSD95_ab'}, ...
    'VGat', {'Synapsin_rs','VGat_ab'}, ...
    'GAD', 'GAD_ab', {'Synapsin_rs', 'GAD_ab'}, ...
    'Gephyrin', {'Synapsin_rs', 'Gephyrin_ab'}, ...
    };

for i = 1:length(vis)
    ds = addvisrow(ds, vis{i});
end

ds.ft = [];
ds.ftname = {};

% integrated brightness features
for c = 1:ds.nch
   cname = ds.chlist{c};
   name = sprintf('%s_ib', cname);
   dat = zeros(ds.ntrain,1);
   for i = 1:ds.ntrain
       chdat = min(ds.sg(i).im(6,5:7,5:7,5:7), ds.sg(i).im(c,5:7,5:7,5:7));
       dat(i) = sum(chdat(:));
   end
   
   ds = addfeature(ds, dat, name);
end

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

ds = addfeature(ds, min([glutdat(:,1) sqrt(sum(glutdat(:,2:3).^2,2))], [], 2), 'Glut_L2prepost');
ds = addfeature(ds, min([gabadat(:,1) sqrt(sum(gabadat(:,2:3).^2,2))], [], 2), 'GABA_L2prepost');

% which training examples to actually plot and test on
ds.trainactive = ds.trainlabelconf > 0.4;

%% Test CV error
testft;

%% Plot a 2D separability scatter plot
figure(1), clf;
feat = {'Glut_L2prepost', 'GABA_L2prepost'};
synplot(ds, feat);
