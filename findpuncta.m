function [coord ib] = findpuncta(ds, channel)

for i = 1:ds.ntrain
    sgim = ds.sg{i}.im;
    se = strel('ball', 5, 5)