% Calculate feature vectors from synaptogram data in ds.sg

ds.trainFeatName = cell(1,1);

% integrated brightness features

ibfeat = zeros(ds.ntrain, ds.nch);
ibfeatname = cell(1,ds.nch);

for c = 1:ds.nch
   cname = ds.chlist{c};
   ibfeatname{c} = sprintf('IB_%s', cname);
   
   for i = 1:ds.ntrain
       chdat = ds.sg(i).im(6,5:7,5:7,5:7) .* ds.sg(i).im(c,5:7,5:7,5:7);
       ibfeat(i,c) = sum(chdat(:));
   end
end

ds.ftname = ibfeatname;
ds.ft = ibfeat;

ds.labelcolors = [0 0 0; 1 0 0; 0 0.5 1; 1 1 0];

figure(1), clf;
colors = ds.labelcolors(ds.trainlabel,:);
gscatter(ibfeat(:,5), ibfeat(:,7), ds.labelnames(ds.trainlabel)', colors);
