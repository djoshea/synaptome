%% Test LDA Cross validation error to evaluate feature utility

% ds.trainactive = ds.trainactive & ds.trainlabel ~= 1 & ds.trainlabel ~= 4;

testname = cell(1);
features = cell(1);
labels = cell(1);

testname{1} = '5-Class LDA';
features{1} = ds.ft(ds.trainactive,:);
labels{1} = ds.trainlabel(ds.trainactive);

testname{2} = '2-Class Glut vs. None';
features{2} = ds.ft(ds.trainactive,:);
labels{2} = (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'glut1'))) | ...
            (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'glut2')));

testname{3} = '2-Class GABA vs. None';
features{3} = ds.ft(ds.trainactive,:);
labels{3} = (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'gaba')));

for i = 1:numel(testname)
    X = features{i};
    Y = labels{i};
    fprintf('Test: %s\n',testname{i});
    trerr = sum(classify(X,X,Y) ~= Y) / length(Y);

    cp = cvpartition(Y, 'k', 3);
    ldaMCR_fn = @(xtrain, ytrain, xtest, ytest) sum(classify(xtest, xtrain, ytrain)~=ytest);
    cvout = crossval(ldaMCR_fn, X, Y, 'partition', cp);
    ldaCVerr = sum(cvout) / sum(cp.TestSize);

    fprintf('\tTraining Set Error: %0.4f\n', trerr);
    fprintf('\tCross Validation Error: %0.4f\n', ldaCVerr);
end