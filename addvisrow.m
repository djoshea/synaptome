function ds = addvisrow(ds, visdat, name)
% adds a row to the visualization structure
% visdat is either a string pointing into the trdlist array of channel
% names or a cell array of three strings to form an RGB list
% name is either specified manually or determined from the visdat structure

if(iscell(visdat))
    for c = length(visdat)+1:3
        % fill empty channels with ''
        visdat{c} = '';
    end
end

ds.vis{end+1} = visdat;

if(~exist('name', 'var')) 
    % determine name from visdat directly
    if(ischar(visdat)) % one channel vis
        name = visdat;
    else
        % 3 channel cell
        name = sprintf('%s/%s/%s', visdat{1}, visdat{2}, visdat{3});
    end 
end

ds.visname{end+1} = char(name);
        