x = randn(100, 1);
y = x + 0.1*randn(100,1);

figure(1), clf;
plot(x,y, 'rx');

[coeff score latent] = princomp([x y]);

figure(2), clf;
plot(score(:,1), score(:,2), 'bx');
