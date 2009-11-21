%% Test LDA Cross validation error to evaluate feature utility

X = ds.ft(ds.trainactive,:);
Y = ds.trainlabel(ds.trainactive);

cp = cvpartition(Y, 'k', 10);
ldaMCR_fn = @(xtrain, ytrain, xtest, ytest) sum(classify(xtest, xtrain, ytrain)~=ytest);
cvout = crossval(ldaMCR_fn, X, Y, 'partition', cp);

ldaCVerr = sum(cvout) / sum(cp.TestSize);

fprintf('Cross Validation Error: %0.4f\n', ldaCVerr);