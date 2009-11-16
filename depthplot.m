% clear ds;

disp('Loading Pivot List...');
ds.fnamepivot = 'gordon300/Synapsin.tif.pivots.xls';
ds.pivots = loadpivotlist(ds.fnamepivot);
ds.npivots = length(ds.pivots.x);

ds.indicatorfnames = {'gordon300/glutAllClassified' , ...
    'gordon300/gabaAllClassified'};

%%
disp('Loading Labels...');
ds.classes = {'glut', 'gaba'};
ds.nc = length(ds.classes);

ds.indic = zeros(ds.npivots, ds.nc);
for c = 1:ds.nc
    fid = fopen(ds.indicatorfnames{c});
    data = textscan(fid, '%d');
    ds.indic(:,c) = data{1};
    fclose(fid);
end

ds.labelnames = {'none', 'glut', 'gaba', 'both'};
ds.nl = length(ds.labelnames);
ds.labelids = 0:ds.nl-1;
ds.label = 1*ds.indic(:,1) + 2*ds.indic(:,2);
ds.showlabels = [2 3];
ds.showlabelcolors = {'r' 'b'};
ds.nlshow = length(ds.showlabels);

%% Cortical Depth Projection

% load info on dimensions from first mask, store for later
ds.maskfnames = {'/Volumes/experiments/Nick/_data/_Gordon 20090630 stacks/WT/Processing/segmentation/YFP RBD2 D2AF unwarp bksub100px Otsu-binary.tif', ...
 '/Volumes/experiments/Nick/_data/_Gordon 20090630 stacks/WT/Processing/segmentation/YFP RBD2 D2AF unwarp bksub100px Otsu-binary expanded-250nm.tif', ...   
'/Volumes/experiments/Nick/_data/_Gordon 20090630 stacks/WT/Processing/segmentation/YFP RBD2 D2AF unwarp bksub100px Otsu-binary expanded-500nm.tif', ...
'/Volumes/experiments/Nick/_data/_Gordon 20090630 stacks/WT/Processing/segmentation/YFP RBD2 D2AF unwarp bksub100px Otsu-binary expanded-1000nm.tif'};
ds.masknames = {'YFP Orig', 'YFP 250 nm', 'YFP 500 nm', 'YFP 1000 nm'};
ds.nm = length(ds.masknames); % number of masks

info = imfinfo(ds.maskfnames{1});
ds.nz = length([info.Height]) - 2;  % skip first and last z-plane
nz = ds.nz;
ds.ny = info(1).Height; ny = ds.ny;
ds.nx = info(2).Width; nx = ds.nx;

% voxel dimensions for density calculations
ds.dz = 0.07; % um
ds.dxy = 0.1; % um

disp(sprintf('Image Dimensions: %.2f x %.2f x %.2f um XYZ', ds.dxy * ds.nx, ds.dxy * ds.ny, ds.dz * ds.nz));

% tilt params
ds.yleft = 788; % offset from top of image to pia at left edge
ds.yright = 512; % offset from top of image to pia at right edge

% radians off of horizontal, positive tilts up towards right
ds.theta = atan( (ds.yleft-ds.yright) / ds.nx );
ds.dmax = (ds.ny - max(ds.yleft, ds.yright)) * cos(ds.theta);

ds.depthfn = @(x,y) double(x)*sin(ds.theta) + (double(y) - ds.yleft)*cos(ds.theta);
ds.depth = ds.depthfn(ds.pivots.x, ds.pivots.y);

% depth bins
ds.dbins = linspace(0, ds.dmax, 50);

% calculate depth bin volume
ds.dbindelta = ds.dbins(2) - ds.dbins(1);
ds.dbindeltay = ds.dbindelta / cos(ds.theta);
ds.dbinvol = (ds.nx/cos(ds.theta)*ds.dbindelta * ds.dxy^2) * (ds.nz*ds.dz);

% calculate depth overlay
[Yvals Xvals] = meshgrid(1:ny, 1:nx);
ds.depthovl = ds.depthfn(Xvals, Yvals); % depth at each point in image


%% By Type Projection

disp('Computing depth histogram...');
cmap = jet(ds.nlshow);

ds.dcounts = zeros(ds.nlshow, length(ds.dbins));
for k = 1:ds.nlshow
    c = ds.showlabels(k);
    ds.dcounts(k,:) = histc(ds.depth(ds.label == ds.labelids(c)), ds.dbins);
end

%% Nick's YFP Masks: for YFP Density calculations and Synaptic Density Profiles within mask
ds.mask = cell(ds.nm,1);

for m = 1:ds.nm
    fprintf('Processing %s Mask...\n', ds.masknames{m});
    ds.mask{m}.dcounts = zeros(ds.nlshow, length(ds.dbins));
    ds.mask{m}.pivotinds = zeros(ds.npivots, 1);
    ds.mask{m}.dpxcounts = zeros(size(ds.dbins));
    
    for z = 2:nz+1 % skip first and last plane (nz already excludes these)
       
        % read in z-th slice of m-th mask
        fprintf('    Reading Slice %d/%d: ', z, nz+1);
        zimg = (imread(ds.maskfnames{m}, z) > 0)'; 
        fprintf('\t%.3f occupied by Mask,', nnz(zimg) / numel(zimg) );
        
        % find all pivots within yfp mask and mark indices in pivotinds
        inds = (ds.pivots.z == z) & (zimg(sub2ind(size(zimg), ds.pivots.x, ds.pivots.y)));
        ds.mask{m}.pivotinds(inds) = ds.mask{m}.pivotinds(inds) + 1; 
        
        fprintf('\t%d/%d Synapses\n', nnz(inds), nnz(ds.pivots.z==z));
        
        % bin mask pixels by depth and add to running tally 
        % (not synapses, just mask volume binned by depth)
        ds.mask{m}.dpxcounts = ds.mask{m}.dpxcounts + histc(ds.depthovl(zimg), ds.dbins)';
        
        if(z == 2) % a quick overlay display
            figure(4), clf;
            imshow(zimg');
            hold on
     
            % plot synapses, colored by whether they are or are not in mask
            plot(ds.pivots.x(ds.mask{m}.pivotinds & ds.pivots.z == z), ds.pivots.y(ds.mask{m}.pivotinds & ds.pivots.z == z), 'r.', 'MarkerSize', 1);
            plot(ds.pivots.x(~ds.mask{m}.pivotinds & ds.pivots.z == z), ds.pivots.y(~ds.mask{m}.pivotinds & ds.pivots.z == z), 'c.', 'MarkerSize', 1);      
%             plot(ds.pivots.x(ds.mask{m} & ds.label == 2 & ds.pivots.z == z), ds.pivots.y(ds.mask{m} & ds.label == 2 & ds.pivots.z == z), 'r.', 'MarkerSize', 2);
%             plot(ds.pivots.x(~ds.mask{m} & ds.label == 2 & ds.pivots.z == z), ds.pivots.y(~ds.mask{m} & ds.label == 2 & ds.pivots.z == z), 'c.', 'MarkerSize', 2);
            
            % plot depth boundaries
            plot([1 ds.nx], [ds.yleft, ds.yright], 'y-', 'LineWidth',2);
            plot([1 ds.nx], [ds.yleft + ds.dmax/cos(ds.theta), ds.yright + ds.dmax/cos(ds.theta)], 'y-', 'LineWidth',2);
            
            hold off
            drawnow
%             pause
        end
    end
    
    % loop over labels and bin masked
    for k = 1:ds.nlshow 
        label = ds.labelids(ds.showlabels(k));
        ds.mask{m}.dcounts(k,:) = hist(ds.depth(ds.label == label & ds.mask{m}.pivotinds > 0), ds.dbins);
    end
end

