%% Create data channels from imaging data
% this involves copying processed data from ds.sg.im to ds.dat
% use addchannel(ds, chdat, chname) to add to ds.trd

% ds.ch is a ndatachannels x ntrain x sgdim 5-D array of processed data
ds.trd = zeros([0 ds.ntrain ds.sgdim]);
ds.trdlist = cell(0);

% for now, just copy raw image data as channels
for i = 1:length(ds.chlist)
    chname = ds.chlist{i};
    ds = addchannel(ds, getimchannel(ds, chname), chname);
end

%% Specify synaptogram visualization structure
% this involves creating a cell array with one element for each row
% each row consists of 1 or 3 channel names that index into ds.ch
% 1 name implies grayscale, 3 names are used as RGB channels

ds.vis = {};
ds.visname = {};

vis = { {'','Synapsin'}, 'Bassoon', ...
    {'Synapsin','VGlut1'}, {'Synapsin', 'VGlut2'}, ...
    {'Synapsin', 'Gephyrin', 'PSD95'}, ...
    'VGat', 'GAD', 'Gephyrin', 'YFP'};

for i = 1:length(vis)
    ds = addvisrow(ds, vis{i});
end

ds.ft = [];
ds.ftname = {};

% integrated brightness features
nch = size(ds.trd,1);
for c = 1:ds.nch
   cname = ds.chlist{c};
   name = sprintf('IB_%s', cname);
   dat = zeros(ds.ntrain,1);
   for i = 1:ds.ntrain
       chdat = min(ds.sg(i).im(6,5:7,5:7,5:7), ds.sg(i).im(c,5:7,5:7,5:7));
       dat(i) = sum(chdat(:));
   end
   
   ds = addfeature(ds, dat, name);
end

figure(1), clf;
feat = {'IB_VGlut1', 'IB_PSD95'};
synplot(ds, feat);
