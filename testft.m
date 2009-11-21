%% Test LDA Cross validation error to evaluate feature utility

% ds.trainactive = ds.trainactive & ds.trainlabel ~= 1 & ds.trainlabel ~= 4;

X = ds.ft(ds.trainactive,:);
Y = ds.trainlabel(ds.trainactive);

trerr = sum(classify(X,X,Y) ~= Y) / length(Y);

cp = cvpartition(Y, 'k', 2);
ldaMCR_fn = @(xtrain, ytrain, xtest, ytest) sum(classify(xtest, xtrain, ytrain)~=ytest);
cvout = crossval(ldaMCR_fn, X, Y, 'partition', cp);
ldaCVerr = sum(cvout) / sum(cp.TestSize);

fprintf('Training Set Error: %0.4f\n', trerr);
fprintf('Cross Validation Error: %0.4f\n', ldaCVerr);