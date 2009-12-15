%% Test LDA Cross validation error to evaluate feature utility

% ds.trainactive = ds.trainactive & ds.trainlabel ~= 1 & ds.trainlabel ~= 4;

testname = cell(1);
features = cell(1);
labels = cell(1);

ft = ds.ft;

testname{1} = '5-Class LDA';
ft = getfeature(ds, {'VGlut1_iab', 'VGlut2_iab', 'Glut_L2pre', 'GABA_L2pre', 'PSD95_iab', 'Gephyrin_iab'});
features{1} = ft(ds.trainactive,:);
labels{1} = ds.trainlabel(ds.trainactive);

testname{2} = '2-Class Glut vs. None';
ft = getfeature(ds, {'Glut_L2pre', 'PSD95_iab'});
features{2} = ft(ds.trainactive,:);
 labels{2} = (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'glut1'))) | ...
             (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'glut2'))) | ...
             (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'glutBoth')));
%labels{2} =  (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'glut')));

testname{3} = '2-Class GABA vs. None';
ft = getfeature(ds, {'GABA_L2pre', 'Gephyrin_iab'});
features{3} = ft(ds.trainactive,:);
labels{3} = (ds.trainlabel(ds.trainactive) == find(strcmp(ds.labelnames, 'gaba')));

errors = zeros(nnz(ds.trainactive), numel(testname));

for i = 1:numel(testname)
    X = features{i};
    Y = labels{i};
    fprintf('Test: %s\n',testname{i});
    errors(:,i) = (classify(X,X,Y) ~= Y);
    trerr = sum(errors(:,i)) / length(Y);

%     cp = cvpartition(Y, 'k', 3);
%     ldaMCR_fn = @(xtrain, ytrain, xtest, ytest) sum(classify(xtest, xtrain, ytrain,'quadratic')~=ytest);
%     cvout = crossval(ldaMCR_fn, X, Y, 'partition', cp);
%     ldaCVerr = sum(cvout) / sum(cp.TestSize);

    fprintf('\tTraining Set Error: %0.4f\n', trerr);
%     fprintf('\tCross Validation Error: %0.4f\n', ldaCVerr);
end

%% Plot a 2D separability scatter plot
figure(1), clf;
feat = {'GABA_L2pre', 'Gephyrin_iab'};
synplot(ds, feat,errors(:,3), {'gaba'});

figure(2), clf;
feat = {'Glut_L2pre', 'PSD95_iab'};
synplot(ds, feat, errors(:,2), {'glut1', 'glut2', 'glutBoth'});
% synplot(ds, feat, errors(:,2), {'glut'});