function [glch, glsrc, gldet, glwl] = map_gs(enum, varargin)
% MAP_GS Construct global spectroscopic enumeration.
%
%   [glch, glsrc, gldet, glwl] = MAP_GS(enum)
%
%   MAP_GS constructs a globally indexed spectroscopic channel table from the provided
%   enumeration, alongside maps from the global indices to the source and detector optode
%   indices
%
%   Paramters:
%
%     enum:   An enumeration returned from lumo_read
%
%   Optional Parameters:
%
%   'group':  An integer specifiying the group index of the enumeration to map. Defaults to
%             group index 1.
%
%   Returns:
%
%     glch:   The global channel list matrix [nc x 3] maps the channel index to the global
%             enumeration.
%
%             glch(ci, 1) -> global source index of channel ci
%             glch(ci, 2) -> global wavelength index of channel ci
%             glch(ci, 3) -> global detector index of channel ci
%
%     glsrc:  The global to local source mapping array [n_src_pos x 1] maps a global source
%             index into the canonical enumeration, allowing e.g., the location of an optode
%             or source power to be determined.
%
%             glsrc(glch(ci,1)).node_idx
%             glsrc(glch(ci,1)).src_idx
%             glsrc(glch(ci,1)).source_optode_idx
%
%     gldet:  The global to local detector mapping array [n_det_pos x 1] maps a global
%             detector index into the canonical enumeration, allowing e.g., the location of 
%             an optode to be looked up in the associated template layout.
%
%             gldet(glch(ci,3)).node_idx
%             glsrc(glch(ci,1)).det_idx
%             gldet(glch(ci,3)).detector_optode_idx
%
%     glwl:   The global wavelength mapping.
%
%             glwl(glch(ci,2)) -> source wavelength of channel ci
%
%   Details:
%
%   The canonical form of the LUMO enumeration represents the system in a node-local format,
%   such that each channel is defined by reference to the local sources and detectors of a
%   given node e.g.,
%
%   channel(i) = ( src_node(j) src_index(k) det_node(l) det_idx(m) )
%
%   Note that in this context, each source is a single emitter at a single wavelength.
%
%   It is common amongst fNIRS and DOT analysis software to represent the system in a global
%   form such that a channel is defined by reference to a set of global source and detector
%   locations, and an emission wavelength, e.g.,
%
%   channel(i) = ( gl_src_index(j) gl_wavelength(k) gl_det_index(l) )
%
%   Whilst less flexible, this organisation is sane insofar as typical instruments (i) have
%   fixed layouts, and (ii) are designed for spectroscopy such that all source locations
%   have emitters of every wavelength.
%
%   The glch matrix provides each of the global indicies for each channel in the
%   enumeration, alongside a global wavelength table indexed by the same. To maintain the
%   association between the global source / detector indices and the physical optodes of the
%   system, two further arrays are provided which link map the indices to the node and
%   optode indices which are indexed in the associated layout.
%
%   For example, if we wish to determine the physical location of the source in channel 4,
%   we might index as follows:
%
%   node_idx = glsrc(glch(2,1));
%   node_id = enum.groups(1).nodes(node_idx).id;
%   src_loc = layout.docks(node_id).coords_3d...;
%
%   Notes:
%
%   Since the globally indexed form has less inherent flexibility than the canonical
%   enumeration used by LUMO, checks are performed to ensure that an appropriate mapping can
%   be performed (e.g., that all source optodes have the same wavelengths, and so on).
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

% First pass over the nodes:
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
n_src_pos = 0;      % Source positions are LUMO optodes of type 'S'
n_det_pos = 0;      % Detector positions are LUMO optodes of type 'D'

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
    
    % These dimesnions are used for preallocatoin of the map
    qimax = max([enum.groups(gi).nodes(ni).srcs.optode_idx]);
    mimax = max([enum.groups(gi).nodes(ni).dets.optode_idx]);
    
    wavelengths = sort(unique([enum.groups(gi).nodes(ni).srcs.wl]));
    
  else
    
    % Assert homogeneity
    assert(length(ql) == ndq);
    assert(length(ml) == ndm);
    assert(qimax == max([enum.groups(gi).nodes(ni).srcs.optode_idx]));
    assert(mimax == max([enum.groups(gi).nodes(ni).dets.optode_idx]));
    assert(all(wavelengths == sort(unique([enum.groups(gi).nodes(ni).srcs.wl]))));
    
  end
  
  % Add the source and detector optode locations
  n_src_pos = n_src_pos + length(ql);
  n_det_pos = n_det_pos + length(ml);
  
end

% Initialise indexing variables
soi = 0;
doi = 0;

% Preallocate map [ max node index x max local optode index ]
%
src_opt_l2g = zeros(n_node, qimax);
det_opt_l2g = zeros(n_node, mimax);

% Second pass over the nodes:
%
% - Loop over all entries in order to build the source and detector to
%   optode mapping.
%
for ni = 1:n_node
  
  % Get the node ID, which forms the link between nodes and docks
  node = enum.groups(gi).nodes(ni);
  nid = node.id;
  
  % Get the source and detector optode indices for this node
  src_optode_idx = unique([enum.groups(gi).nodes(ni).srcs.optode_idx]);
  det_optode_idx = unique([enum.groups(gi).nodes(ni).dets.optode_idx]);
  
  % Add source entries
  for qi = 1:length(src_optode_idx)
    soi = soi + 1;
    assert(src_opt_l2g(ni, src_optode_idx(qi)) == 0);
    
    % Source optode local [nid, optode_idx] -> [gl src opt index]
    src_opt_l2g(ni, src_optode_idx(qi)) = soi;
    src_opt_g2l(soi) = struct('node_idx', ni, 'optode_idx', src_optode_idx(qi));
  end
  
  % Add detector entries
  for mi = 1:length(det_optode_idx)
    doi = doi + 1;
    assert(det_opt_l2g(ni, det_optode_idx(mi)) == 0);
    
    % Detector optode local [nid, optode_idx] -> [gl det opt 1-index]
    det_opt_l2g(ni, det_optode_idx(mi)) = doi;
    det_opt_g2l(doi) = struct('node_idx', ni, 'optode_idx', det_optode_idx(mi));
  end
  
end

% Build MeasList
%
% Loop over each channel, forming the global mapping based upon the
nc = length(enum.groups(gi).channels);
glch = zeros(nc, 3);

for ci = 1:nc
  
  ch = enum.groups(gi).channels(ci);             % The channel
  ndq = enum.groups(gi).nodes(ch.src_node_idx);  % The source node
  ndm = enum.groups(gi).nodes(ch.det_node_idx);  % The detector node
  
  ml_qi = src_opt_l2g(ch.src_node_idx, ndq.srcs(ch.src_idx).optode_idx);
  ml_mi = det_opt_l2g(ch.det_node_idx, ndm.dets(ch.det_idx).optode_idx);
  ml_wi = find(ndq.srcs(ch.src_idx).wl == wavelengths);
  
  glch(ci,1) = ml_qi;
  glch(ci,2) = ml_wi;
  glch(ci,3) = ml_mi;
  
end

glsrc = src_opt_g2l;
gldet = det_opt_g2l;
glwl = wavelengths;

end

