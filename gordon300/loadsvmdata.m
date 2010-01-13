% For Gordon200 set Based off GoldSVM1n2 file

% function dat = loadsvmdata(fname)
fname = 'gordon200\GoldSVM1n2';
fid = fopen(fname);
fseek(fid, 0, 'bof');

% find dimension of feature space
fscanf(fid, '%d', 1);
ln = fgetl(fid); % get remainder of line after first class label
k = textscan(ln, '%d:%f');
d = length(k{1});
fseek(fid, 0, 'bof');

% read in data
fmatstr = ['%d ' repmat('%*d:%f ', 1, d)];
data = textscan(fid, fmatstr);
fclose(fid);

Y = data{1};
X = [data{2:end}];

syns = Y > 0;

Xsyn = X(syns,:);
Ysyn = Y(syns);

% spherize Xsyn
for j = 1:size(X,2);
    X(:,j) = X(:,j) - mean(X(:,j));
    X(:,j) = X(:,j) / norm(X(:,j));
end

for j = 1:size(Xsyn,2);
    Xsyn(:,j) = Xsyn(:,j) - mean(Xsyn(:,j));
    Xsyn(:,j) = Xsyn(:,j) / norm(Xsyn(:,j));
end

nons = Ysyn == 0;
excs = Ysyn == 1;
inhs = Ysyn == 2;

d = size(Xsyn, 2);
m = size(Xsyn, 1);

% PCA Projection

[coeff score latent] = princomp(X);

figure(1), clf;
% gscatter(score(:,1), score(:,2), Ysyn, 'br');
plot3(score(excs,1), score(excs,2), score(excs, 3), 'g.', 'MarkerSize', 10);
hold on
plot3(score(inhs,1), score(inhs,2), score(inhs, 3), 'r.', 'MarkerSize', 10);
% plot3(score(nons,1), score(nons,2), score(nons, 3), 'b.', 'MarkerSize', 10);

xlabel('PC 1')
ylabel('PC 2');
zlabel('PC 3');
view([0 90]);

box on
grid on

%% LDA Projection

[Xrot rot] = ldaproj(X, Y);

figure(2), clf;
gscatter(Xrot(:,1), Xrot(:,2), Y, 'rgby', '.', 10);
% gscatter(score(:,1), score(:,2), Ysyn, 'br');
% plot3(Xrot(excs,1), Xrot(excs,2), Xrot(excs,3), 'g.', 'MarkerSize', 10);
% hold on
% plot3(Xrot(inhs,1), Xrot(inhs,2), Xrot(inhs,3), 'r.', 'MarkerSize', 10);
% % plot3(score(nons,1), score(nons,2), score(nons, 3), 'b.', 'MarkerSize', 10);
% view(0,90);
% xlabel('PC 1')
% ylabel('PC 2');
% zlabel('PC 3');
title('Fisher Projection');

% box on
% grid on

% %% Quadratic feature space expansion
% 
% Xsynq = zeros(m, d^2 + d);
% for i = 1:length(Ysyn)
%     xr = [1 Xsyn(i,:)];
%     xq = xr'*xr;
%     Xsynq(i,:) = reshape(xq(:,2:end), 1, []);
% end
% 
% [coeff score latent] = princomp(Xsynq);
% 
% figure(3), clf;
% % gscatter(score(:,1), score(:,2), Ysyn, 'br');
% plot3(score(excs,4), score(excs,2), score(excs, 3), 'g.', 'MarkerSize', 10);
% hold on
% plot3(score(inhs,4), score(inhs,2), score(inhs, 3), 'r.', 'MarkerSize', 10);
% plot3(score(nons,4), score(nons,2), score(nons, 3), 'b.', 'MarkerSize', 10);
% 
% xlabel('PC 1')
% ylabel('PC 2');
% zlabel('PC 3');
% grid on
% box off
% 
% %% Expanded feature space LDA projection
% 
% priors = arrayfun(@(k) nnz(Ysyn==k) / m, [1 2]);
% 
% meanfn = @(k) sum(Xsynq(Ysyn==k,:)) / nnz(Ysyn==k);
% means = [meanfn(1); meanfn(2)];
% 
% sigma = zeros(size(Xsynq,2));
% for i = 1:m
%     x = Xsynq(i,:);
%     k = Ysyn(i);
%     sigma = sigma + (x-means(k,:))'*(x-means(k,:));
% end
% sigma = sigma / (m-2);
% 
% M = means;
% W = sigma;
% 
% Mstar = M*W^(-1/2);
% 
% Bstar = cov(Mstar);
% 
% [V D] = eig(Bstar);
% 
% %%
% Xrotq = (V(1:3,:)*Xsynq')';
% 
% %%
% figure(4), clf;
% % gscatter(score(:,1), score(:,2), Ysyn, 'br');
% plot3(Xrotq(excs,1), Xrotq(excs,2), Xrotq(excs,3), 'g.', 'MarkerSize', 10);
% hold on
% plot3(Xrotq(inhs,1), Xrotq(inhs,2), Xrotq(inhs,3), 'r.', 'MarkerSize', 10);
% % plot3(score(nons,1), score(nons,2), score(nons, 3), 'b.', 'MarkerSize', 10);
% view(0,90);
% xlabel('PC 1')
% ylabel('PC 2');
% zlabel('PC 3');
% 
% title('Expanded Feature Space LDA Projection');
% 
% box on
% grid on

