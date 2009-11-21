function [feat idx] = getfeature(ds, names)

if(ischar(names))
    names = {names};
end
N = length(names);
idx = zeros(1,N);
feat = zeros(ds.ntrain, N);

for i = 1:length(names)
    idxs = find(strcmp(names{i}, ds.ftname));
    if(~isempty(idxs))
        idx(i) = idxs(1);
        feat(:,i) = ds.ft(:, idx(i));
    else
        error(sprintf('Could not find feature "%s"', names{i}));
    end
end


