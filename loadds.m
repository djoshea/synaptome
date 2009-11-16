dsname = 'gordon300';

% All synapsin pivot locations
fnamePivots = strcat(dsname, '/Synapsin.tif.pivots.xls');

% Training synapsin pivot locations
fnameTrainPivots = strcat(dsname, '/Synapsin.tifSubsetObjects.xls');

% Training synapsin pivot human labels
fnameHumanLabels = strcat(dsname, '/classifications.mat');

% Indicator Labels from SVM classifier for full dataset
indicatorfnames = {'gordon300/glutAllClassified' , ...
    'gordon300/gabaAllClassified'};

% Output classes
ds.classes = {'glut', 'gaba'};

% Labels for combinatorial patterns
ds.labelnames = {'none', 'glut', 'gaba', 'both'};
ds.labelcolors = {'k', 'r', 'c', 'm'};
ds.showlabels = [2 3]; % labels to show
ds.showlabelcolors = {'r' 'c'}; % with these colors

% Channel List for Synaptogram
ds.chlist = {'Bassoon', 'DAPI', 'GAD', 'Gephyrin', 'PSD95', 'Synapsin', 'VGat', 'VGlut1', 'VGlut2', 'YFP'};
ds.nch = length(ds.chlist);
ds.sgdim = [11 11 11]; % dimensions of synaptogram

%% Train Pivot List
disp('Loading Training Pivot List...');
ds.trainPivots = loadpivotlist(fnameTrainPivots);
ds.ntrain = length(ds.trainPivots.x);

disp('Loading Human Labels...');
load(fnameHumanLabels);
ds.votes = classifications;

ds.nvotes = size(ds.votes,3);
ds.nvotecols = size(ds.votes, 2);
% replace two column binary voting with ordinal number that indexes into
% labels
mult = repmat(arrayfun(@(x) 2^(x-1), 1:ds.nvotecols), [ds.ntrain, 1, ds.nvotes]);
ds.votesnumeric = squeeze(sum(ds.votes .* mult, 2) + 1); % ntrain x nvotes
ds.trainlabel = mode(ds.votesnumeric, 2); % most frequent vote, tie goes to smaller index
ds.trainlabelrunnerup = zeros(ds.ntrain, 1);
ds.trainlabelconf = zeros(ds.ntrain,1);
for i = 1:ds.ntrain
    row = ds.votesnumeric(i,:);
    ds.trainlabelrunnerup(i) = mode(row(row~=ds.trainlabel(i)));
    ds.trainlabelconf(i) = (nnz(row == ds.trainlabel(i)) - ...
        nnz(row == ds.trainlabelrunnerup(i))) * 1/ds.nvotes;
end

disp('Loading Training Synaptogram Data...');
ds.sg = cell(ds.ntrain, 1);
for i = 1:ds.ntrain
    disp(sprintf('    Synaptogram %4d', i));
    ds.sg{i} = loadsg(ds,i); 
end

% convert to struct array for easy access
st = struct([]);
for i = 1:ds.ntrain
    st(i).i = ds.sg{i}.i;
    st(i).x = ds.sg{i}.x;
    st(i).y = ds.sg{i}.y;
    st(i).z = ds.sg{i}.z;
    st(i).im = ds.sg{i}.im;
    st(i).coordstr = ds.sg{i}.coordstr;
end
ds.sg = st;