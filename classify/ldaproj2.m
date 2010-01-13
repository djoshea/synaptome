function [w xval] = ldaproj2(X,Y) 

m = size(X,1);
d = size(X,2);
ks = unique(Y);
K = length(ks);

allmean = mean(X);

means = zeros(K, d);
for k = 1:K
    means(k,:) = sum(X(Y==ks(k),:)) / nnz(Y==ks(k));
end

Sb = zeros(d);
for c = 1:K
    Sb = Sb + (means(c,:) - allmean)'*(means(c,:) - allmean);
end

Sw = zeros(d);
for i = 1:m
    x = X(i,:);
    mc = means(find(ks==Y(i)),:);
    Sw = Sw + (x-mc)'*(x-mc);
end

[U S V] = svd(Sb);
sbhalf = U*sqrt(S)*V';

[V D] = eig(sbhalf*Sw^(-1)*sbhalf);

w = sbhalf^(-1) * V(:,1);

xval = X*w;
