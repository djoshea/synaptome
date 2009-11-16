Y = (-150:149 >= 0)';
m = 300;
Xnoise = 500*randn(m, 5);
Xsig = [50*(2*Y-1) -5*(2*Y-1)] + 1*randn(m,2);
X = [Xsig Xnoise];

randrot = 5*randn(size(X,2));
% randrot = eye(size(X,2));
X = (randrot*X')';

% spherize Xsyn
for j = 1:size(X,2);
    X(:,j) = X(:,j) - mean(X(:,j));
    X(:,j) = X(:,j) / norm(X(:,j));
end

% Original Space
figure(1), clf;
gscatter(X(:,1), X(:,2), Y, 'rgby', '.', 10);
title('Original Feature Coordinates');

% PCA projection

[score coeff latent] = princomp(X);

figure(2), clf;
gscatter(coeff(:,1), coeff(:,2), Y, 'rgby', '.', 10);
title('PCA Projection');

% LDA Projection
[Xrot rot] = ldaproj(X,Y);

figure(3), clf;
gscatter(Xrot(:,1), Xrot(:,2), Y, 'rgby', '.', 10);
title('Fisher Projection');