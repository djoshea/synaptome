Y = 2*(-150:149 >= 0)' - 1;
m = length(Y);

X = zeros(length(Y), 3);
X(:,1) = randn(length(Y), 1);
X(:,2) = 0.1*randn(length(Y),1) + X(:,1) + Y;
X(:,3) = 0.1*randn(length(Y),1);

% spherize X
% for j = 1:size(X,2);
%     X(:,j) = X(:,j) - mean(X(:,j));
%     X(:,j) = X(:,j) / norm(X(:,j));
% end

% Original Space
figure(1), clf;
gscatter(X(:,1), X(:,2), Y, 'rgby', '.', 10);
title('Original Feature Coordinates');

bins = linspace(-0.3, 0.3, 40);

% PCA projection

[score coeff latent] = princomp(X);

% spherize X
for j = 1:size(coeff,2);
    coeff(:,j) = coeff(:,j) - mean(coeff(:,j));
    coeff(:,j) = X(:,j) / norm(coeff(:,j));
end

cpcaP = histc(coeff(Y==1,1), bins);
cpcaN = histc(coeff(Y==-1,1), bins);

figure(2), clf;
plot(bins, cpcaP, 'b-', 'LineWidth', 2);
hold on
plot(bins, cpcaN, 'r-', 'LineWidth', 2);
title('PCA Projection');

% LDA Projection
[w xval] = ldaproj2(X,Y);
xval = xval - mean(xval);
xval = xval / norm(xval);

cldaP = histc(xval(Y==1,1), bins);
cldaN = histc(xval(Y==-1,1), bins);

figure(3), clf;
plot(bins, cldaP, 'b-', 'LineWidth', 2);
hold on
plot(bins, cldaN, 'r-', 'LineWidth', 2);
title('LDA Projection');
