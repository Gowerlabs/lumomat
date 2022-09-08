function [layout] = read_layout(fn, varargin)
% LUMOFILE.READ_LAYUOT Read a LUMO layout file
%
% [layout] = LUMOFILE.READ_LAYOUT(filename) 
%
% LUMOFILE.READ_LAYOUT reads a standard LUMO layout file in JSON format for disc, returning
% an output structure in standard format.
%
%   Paramters:
%
%   'filename': The path of the LUMO layout file to load.
%                           
%   Optional Parameters:
%
%   Returns:
%
%     layout:   An structure containing the following fields
%
%               id:        a unique identifier
%               name:      friendly name for the group
%               docks:     an array of dock description structures
%               dockmap:   a mapping from dock id to the docks array indices
%               landmarks: (optional) an array of structures of named points on the laytout
%               dims_2d:   (optional) the bounding dimensions of the 2D layout
%               dims_3d:   (optoinal) the bounding dimensions of the 3D layout
%
%   Details:
%
%   The role of the layout file in a LUMO system is to provide a physical reference for each
%   of the channels which are formed in a given group of tiles. The layout file enumerates
%   each dock present in a given group (which could physically be a single cap, headband, or
%   patch).
%
%   Docks:
%
%   Each dock has a numeric identifier, e.g., layout.docks(i).id, which is unique within
%   the layout. When a LUMO tile is present in a dock, it assumes this ID, and thus this
%   forms the link between the logical channels of the system and the physical locations.
%
%   To enable easy lookup of the dock array from the node ID, a mapping array 'dockmap' is
%   provided which maps from ID to array index.
%
%   In addition to the ID, each dock contains an array of optodes descriptors. Each optode
%   is defined by a structure:
%
%     optode.name:       A name for the optode (not used for indexing)
%     optode.coords_2d:  Co-ordinates of the optode in a flattened 2D representation
%     optode.coords_3d:  Co-orindates of the optode in 3D space
%
%   By navigating these structures in conjunction with the canonical channel listing built
%   with LUMOFILE.READ_LUMO, one may construct a complete physical representaiton of the
%   system.
%
%   Landmarks:
%
%   When the LUMO system is to be used as part of DOT reconstruction procedure, it is often
%   neccesary to register the location of the optodes with reference to physiological
%   landmarks. Landmarks are defined by a name (e.g. 'inion') and the physical location in
%   the same 3D co-ordinate space as the optodes.
%
%
%   (C) Gowerlabs Ltd., 2022
%

if exist(fn, 'file') ~= 2
  error('Layour file (%s) cannot be found\n', fn)
end

try
  layout_raw = jsondecode(fileread(fn));
catch e
  fprintf('Error parsing layout file %s\n', fn);
  rethrow(e);
end

layout = lumofile.proc_layout(layout_raw);
