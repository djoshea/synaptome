function ds = addfeature(ds, ftdat, ftname)

ds.ft(:,end+1) = ftdat;
ds.ftname{end+1} = char(ftname);