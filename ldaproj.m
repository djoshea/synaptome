function [Xrot rot] = ldaproj(X,Y) 

m = size(X,1);
d = size(X,2);
ks = unique(Y);
K = length(ks);

priors = arrayfun(@(k) nnz(Y==k) / m, 1:K);

means = zeros(K, d);
for k = 1:K
    means(k,:) = sum(X(Y==ks(k),:)) / nnz(Y==ks(k));
end

sigma = zeros(d,d);
for i = 1:m
    x = X(i,:);
    kind = find(ks==Y(i));
    sigma = sigma + (x-means(kind,:))'*(x-means(kind,:));
end
sigma = sigma / (m-2);

M = means;
W = sigma;
Mstar = M*W^(-1/2);
Bstar = cov(Mstar);
[V D] = eig(Bstar);

Vl = W^(-1/2)*V';
rot = Vl';

Xrot = (rot*X')';
