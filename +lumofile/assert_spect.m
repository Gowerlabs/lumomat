function [is_spect] = assert_spect(enum, varargin)
% ASSERT_SPECT Assert that an enumeration confirms to a spectroscopic scheme
%
%   [is_spect] = ASSERT_SPECT(enum)
%
%   ASSERT_SPECT asserts that an enumeration can be represented in a spectroscopic scheme
%   whereby every source optode contains the same illumination wavelengths.
%
%   Paramters:
%
%     enum:     An enumeration returned from lumo_read
%
%   Returns:
%
%     is_spect: Boolean indicating if the enumeration conforms to requirements
%
%
%   (C) Gowerlabs Ltd., 2022
%

% Parse options
p = inputParser;
addOptional(p, 'group', 1, @(x) (isnumeric(x) && x > 0));
parse(p, varargin{:});

% Set group ID
gi = p.Results.group;

% Check the group index is sensible
if (gi > length(enum.groups))
  error('Requested group %d exceeds that available in the enumeration',gi)
end

% Set true and change if we fail an assertion
is_spect = true;

% Iterate over each node
%
% - Determine total number of 'source locations' and 'detector locations' in NIRS parlance.
%   The underlying assumption is that two tiles do not share a source or detector optode,
%   which is always true.
%
% - Enumerate all wavelengths in the system
%
% - Establish that the enumeration describes an homogenous set of nodes each with the same
%   configuration of sources, detectors, and optodes, in order that the mapping can be
%   constructed.
%

n_node = length([enum.groups(gi).nodes]);

for ni = 1:n_node
  
  % Get a list of all source and detector optode indices
  ql = unique([enum.groups(gi).nodes(ni).srcs.optode_idx]);
  ml = unique([enum.groups(gi).nodes(ni).dets.optode_idx]);
  
  if ni == 1
    
    % On the first iteration, assign the variables in order that we can
    % subsequent assert homogeneity.
    ndq = length(ql);
    ndm = length(ml);  
    wavelengths = [enum.groups(gi).nodes(ni).srcs.wl];    
    
  else
    
    % Assert homogeneity
    try
      assert(length(ql) == ndq);
      assert(length(ml) == ndm);
      assert(all(wavelengths == [enum.groups(gi).nodes(ni).srcs.wl]));
    catch e
      is_spect = false;
    end
      
  end
    
end