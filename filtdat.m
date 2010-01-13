function [filt info] = filtdat(ds, dat, params)
% returns a filtered version of data
% dat is either a channel name (in trdlist) or a ntrain x sgdim... array
% params is a struct with field type and other type-specific settings
% filt is an array ntrain x sgdim1 x ...

% Performs some kind of filtering or computational operation on some
% 	imaging data and returns the results and optionally some metadata.
% 	Filtdat merely returns the filtered output and metadata, it is not 
%   designed to make any changes directly to ds. This enables filters to be
%   chained easily before adding the output to the trd list via addchannel.
%   The filtered output may be easily added back to ds.trd via:
%
%   filtparams = struct('type', 'typename', 'option1', value1);
%   ds = addchannel(ds, filtdat(ds, 'InputChannelName', filtparams), 'OutputChannelName');
% 
% 	dat is the primary input to the operation, specified either directly
% 	  as an array ntrain x sgdim or by name (to be retrieved using
% 	  getchannel)
% 	params is a struct specifying the details of the operation to perform
% 	  params.type specifies which operation to perform and must be 
% 	  specified. Other values under params depend on the operation being
% 	  performed. These values can include parameters used for all synapses
% 	  globally, arrays of parameters used for individual synapses, or even
% 	  entire additional imaging data arrays.
% 	filt is the primary output of the operation as an array ntrain x sgdim
% 	info is a structure containing metadata returned by some operations.
% 	  These fields may contain arrays of values generated for each synapse
% 	  individually or secondary output channel arrays.
% 	 
%    Invidual filters are typically implemented as typename_filt.m files in
%    the filters directory (but may be placed anywhere on the matlab path).
%    They may also be implemented within this file where indicated.
%
% 	 See the documentation for each invidual filter function for more specific
% 	 guidance on usage.

info = [];
filt = [];

if(~exist('params', 'var') || ~isfield(params, 'type'))
    error('Filter type not specified in params.type');
end

if(ischar(dat))
    % get data channel by name
    orig = getchannel(ds, dat);
else
    % assume dat is the data
    orig = dat;
end

% If we find the corresponding _filt.m file on the path, run it
mfilename = sprintf('%s_filt.m', params.type);
funcname = sprintf('%s_filt', params.type);
if(exist(mfilename, 'file'))
    eval(sprintf('[filt info] = %s(ds, orig, params);', funcname));
    return;
end

% can choose to implement other filters here, or create a file called
% type_filt.m on the matlab path or in the filters directory where type ==
% params.type. if several filters are implemented, be sure to end each
% filter's if(params.type == 'mytype') block with return to avoid the
% terminal error message


% if we haven't returned by now, the filter does not exist
error(sprintf('Filter "%s" not implemented', params.type));

