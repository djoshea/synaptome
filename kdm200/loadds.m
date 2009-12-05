ds = [];
ds.dsname = 'kdm200';

% Training synapsin pivot human labels
fnameHumanLabels = strcat('humanlabels.mat');

% Output classes
ds.classes = {'glut1', 'glut2', 'gaba'}; % voting column names

% Labels for combinatorial patterns
ds.labelnames = {'none', 'glut1', 'glut2', 'glutBoth', 'gaba'};
ds.labelmatch = [0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1]; % vote patterns each label matches
ds.labelcolors = [0 0 0; 1 0 0; 0.5 0 0; 1 1 0; 0 0.5 1]; % color for each label
ds.nlabels = numel(ds.labelnames);

% Channel List for Synaptogram
ds.chlist = {'Synapsin', 'Bassoon', 'VGlut1', 'VGlut2', 'PSD95', 'GAD', 'VGat', 'Gephyrin'};
ds.nimch = numel(ds.chlist);
ds.sgdim = [11 11 11]; % dimensions of synaptogram (Z Y X)
ds.sgaspect = [70 100 100]; % aspect ratio of synaptogram (Z Y X)

disp('Loading Human Labels...');
load(fnameHumanLabels);
ds.votes = votes;
ds.nvotes = size(ds.votes,3);
ds.nvotecols = size(ds.votes, 2);
ds.ntrain = size(ds.votes,1);

% replace columnwise binary voting with indices into ds.labelnames array
% based on ds.labelmatch
ds.votesbylabel = zeros(ds.ntrain, ds.nlabels);
for i = 1:ds.ntrain % loop over pivots
    for v = 1:ds.nvotes % loop over voters
        for n = 1:ds.nlabels % loop over label templates
            ds.votesbylabel(i,n) = ds.votesbylabel(i,n) + ...
                all(votes(i,:,v) == ds.labelmatch(n,:));
        end
    end
end

% compute winner and runner-up, and from the difference a confidence metric
[val ds.trainlabel] = max(ds.votesbylabel, [], 2); % most frequent vote, tie goes to smaller index
ds.trainlabelrunnerup = zeros(ds.ntrain, 1);
ds.trainlabelconf = zeros(ds.ntrain,1);
for i = 1:ds.ntrain
    row = ds.votesbylabel(i,:);
    row(ds.trainlabel(i)) = 0;
    [val ds.trainlabelrunnerup(i)] = max(row,[],2);
    if(ds.trainlabelrunnerup(i) == ds.trainlabel(i))
        % received all the votes and happened to be first
        ds.trainlabelconf(i) = 1;
    else
        ds.trainlabelconf(i) = (ds.votesbylabel(i,ds.trainlabel(i)) - ...
            ds.votesbylabel(i,ds.trainlabelrunnerup(i))) / ds.nvotes;
    end
end

% adjustments to votes
ds.trainlabel(168) = 4;
ds.trainlabel(81) = 4;
ds.trainlabel(68) = 3;
ds.trainlabel(71) = 3;
ds.trainlabel(62) = 2;
ds.trainlabel(120) = 5;
ds.trainlabel(28) = 1;


disp('Loading Training Synaptogram Data...');
sgdir = 'U:\Brad\KDM 090416b\Decision tree experiment\puncta';

ds.sg = cell(ds.ntrain, 1);
chlistabbrev = {'01syn', '02bas', '03vglut1', '04vglut2', '05psd', ...
    '06gad', '07vgat', '08geph'};

for i = 1:ds.ntrain
    disp(sprintf('    Synaptogram %4d', i));
    sg = [];
    sg.im = zeros(ds.nimch, ds.sgdim(3), ds.sgdim(2), ds.sgdim(1)); % ch, z, y, x
    sg.str = ds.labelnames{ds.trainlabel(i)};

    for c = 1:ds.nimch
        fname = sprintf('%s\\%s\\T%05d.tif', sgdir, chlistabbrev{c}, i);
        for z = 1:ds.sgdim(3)
            sg.im(c,z,:,:) = imread(fname, z);
        end
        
        % for some reason, central pixel at 6, 6, 6 is set to 255
        % reset it to the 3x3 neighborhood average (excluding me) for now
        nhood = squeeze(sg.im(c,5:7,5:7,5:7));
        nhood(2,2,2) = Inf;
        nhood = nhood(nhood ~= Inf);
        sg.im(c,6,6,6) = mean(nhood);
    end
    ds.sg{i} = sg;
end

sg = ds.sg;

% convert to struct array for easy access
st = struct([]);
for i = 1:ds.ntrain
    st(i).i = i;
    st(i).im = ds.sg{i}.im;
    st(i).str = ds.sg{i}.str;
end
ds.sg = st;

