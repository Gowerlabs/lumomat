function [enum, data, events] = read_lumo(lf_dir, varargin)
% LUMOFILE.READ_LUMO Read a LUMO file from disk
%
% [enum, data, events] = LUMOFILE.READ_LUMO(filename)
%
% LUMOFILE.READ_LUMO reads a LUMO file from disk, returning a set of data structures which
% contain a complete description of the system, and the available data.
%
%   Paramters:
%
%   'filename':       The path of the LUMO file to load.
%
%   Optional Parameters:
%
%   'layout':         When a layout is specified, the embedded layout in the specified LUMO
%                     file, if present, is ignored, and the alternative specification
%                     applied. An alternative layout can be specified as either:
%
%                       string: the filename of a valid lumo layout file, in JSON format
%                       struct: a layout structure in the format returned by
%                               lumofile.read_layout (see function help for details)
%
%                     If the provided layout has been constructed by the user, for example
%                     based upon measurements of a physical layout, entries must be present
%                     for each occupied dock in the recording with the apprporiate dock ID.
%
%   Returns:
%
%     enum:   An enumeration of the system containing:
%
%             enum.hub:     a description of the LUMO Hub used for recording
%             enum.groups:  an array of structures describing each group (cap) connected.
%                           LUMO files only store a single group, so this array is of length 1
%                           and indexing is not required. The form of this structure is
%                           described further below.
%
%     data:   A structure of data form the selected group.
%
%     events: An array of strcutures details events recorded during recording.
%
%   The enumeration (enum)
%
%   The canonical representation of a LUMO system uses a node local indexing format for its
%   channel descriptors. For example, a channel can be defined as the tuple:
%
%   (src_node_idx(j), src_idx(k), det_node_idx(l), det_idx(m))
%
%   This informaiton is exposed in the retuned enumeration:
%
%   >> ch = enum.groups(gi).channels(98)
%
%   ch =
%
%   struct with fields:
%
%     src_node_idx: 1
%          src_idx: 2
%     det_node_idx: 1
%          det_idx: 2
%
%   One may inspect the nature of the sources or detectors with reference to the node to
%   which it belongs, e.g., the wavelength of a source index 2 on node index 1:
%
%   nodes = enum.groups.nodes;
%   wl = nodes(ch.src_node_idx).srcs(ch.src_idx).wl
%
%   wl =
%
%      735
%
%   The physical location of a source is determined by its association with an optode:
%
%   >> optode_idx = nodes(ch.src_node_idx).srcs(ch.src_idx).optode_idx
%
%   optode_idx =
%
%   6
%
%   And the actual physical location of an optode is determined by an associated layout
%   structure. The docks of a layout are linked to enumeration by the node ID. For
%   convenience, the layout structure contains a map from the ID to the index:
%
%   layout = enum.groups.layout;
%   node_id = nodes(ch.src_node_idx).id;
%
%   >> optode = layout.docks(layout.dockmap(node_id)).optodes(optode_idx)
%
%   optode =
%
%     struct with fields:
%
%         name: 'B'
%     coords_2d: [1×1 struct]
%     coords_3d: [1×1 struct]
%
%
%   Note that a 'node' is synonymous with a LUMO tile, or a module in the nomencalature of
%   other fNIRS/DOT formats.
%
%   The canonical representation of a LUMO system allows for a complete representation of
%   an arbitrary LUMO system.
%
%   The data:
%
%   The primary purpose of the data structure is to provide the channel intensity data and
%   fixed associated metadata.
%
%   The data are returned as an array of structures for each group referred to in the
%   enumeration, though since current files only contain a single group, indexing of the
%   array can be dropped.
%
%   data.nchns:     (nc) the number of channels in the recording of this group
%   data.nframes:   (nf) the number of frames in the recording of this group
%   data.chn_fps:   the number of frames per second
%   data.chn_dt:    the length of a frame in milliseconds
%   data.chn_dat:   [nc x nt] matrix of nc channel intensity measurements over nf frames
%
%
%   (C) Gowerlabs Ltd., 2022
%

%%% TODOS
%
% - Add optode filtering in the style of lufr load
%

% Normalise strings
lf_dir = convertStringsToChars(lf_dir);

ts_load = tic;

% Parse inputs
p = inputParser;
addParameter(p, 'layout', []);
parse(p, varargin{:});

layout_override = p.Results.layout;

% Fix the group for LUMO files
gi = 1;

% Load the file description
lf_desc = load_lumo_desc(lf_dir);

% Construct the canonical enumeration and get the data parameters
%
% Note: this function returns the data parameters in addition to the enumeration in order
%       that the relevant metadata files are only parsed once. It may be prudent to split
%       this function whenever the metadata format is revised to respect this distintion.
fprintf('Constructing canonical enumeration...\n');
[enum, lf_dataparam] = load_lumo_enum(lf_dir, lf_desc);

% Load event markers
if lf_desc.has_events
  events = load_lumo_events(lf_dir, lf_desc);
else
  events = [];
end

% Load layout or take from the user
if isempty(layout_override)
  
  % The user has nor provided a layout, we must use the embedded data if it exists
  if lf_desc.has_layout
    
    try
      enum.groups(gi).layout = lumofile.read_layout(fullfile(lf_dir, lf_desc.lo_fn));
    catch e
      fprintf('LUMO file (%s) invalid: error parsing layout file %s\n', lf_dir, lf_desc.lo_fn);
      rethrow(e)
    end
    
  else
    enum.groups(gi).layout = [];
    
    install_layout_path = lumofile.find_install_layout(enum.groups(gi).id);
    [lnhex, lnname] = lumomat.norm_gid(enum.groups(gi).id);
        
    if isempty(install_layout_path)
    
      warning([...
        'The specified LUMO file (group %s / %i) does not contain an embedded layout file, and no layout '...
        'has been specified when calling this function. The returned enumeration will '...
        'lack layout information, and it will not be possible to convert this file to '...
        'formats which require a layout. Specify an appropriate layout file to suppress '...
        'this warning.'], ...
        lnname, str2num(lnhex));

    else
       warning([...
        'The specified LUMO file (group %s / %i) does not contain an embedded layout file, and no layout '...
        'has been specified when calling this function. The returned enumeration will '...
        'lack layout information, and it will not be possible to convert this file to '...
        'formats which require a layout. Specify an appropriate layout file to suppress '...
        'this warning.\n\n'...
        'NOTE: The layout file for this group has been found at: %s.'], ...
        lnname, str2num(lnhex), install_layout_path);
      
    end
      
      
  end
  
else
  
  % The user has supplied a layout file
  if ischar(layout_override)
    
    try
      enum.groups(gi).layout = lumofile.read_layout(layout_override);
      fprintf('LUMO file using user-specified layout file\n');
    catch e
      fprintf(2, 'An error occurred loading the specified layout file %s\n', layout_override);
      rethrow e
    end
    
  elseif isstruct(layout_override)
    enum.groups(gi).layout = layout_override;
    fprintf('LUMO file using user-specified layout structure (not validated)\n');
  else
    error('The specified layout is neither a layout filename nor structure, consult help');
  end
  
end

if ~isempty(enum.groups(gi).layout)
  fprintf('LUMO file assigned layout (group %d) contains %d docks, %d optodes\n', ...
    gi, length(enum.groups(gi).layout.docks), length(enum.groups(gi).layout.docks)*7);
end

[chn_dat] = load_lumo_data(lf_dir, lf_desc, lf_dataparam);

% LUMO file versions 0.4.0&0.5.0 don't have channel saturation data
if all(lf_desc.lfver == [0 4 0]) || all(lf_desc.lfver == [0 5 0])
    data = struct('chn_dat', chn_dat, ...
      'chn_fps', lf_dataparam.chn_fps, ...
      'chn_dt',  round((1/lf_dataparam.chn_fps)*1000), ...
      'nframes', lf_dataparam.nframes, ...
      'nchns',   lf_dataparam.nchns);
else
    data = struct('chn_dat', chn_dat, ...
          'chn_fps', lf_dataparam.chn_fps, ...
          'chn_dt',  round((1/lf_dataparam.chn_fps)*1000), ...
          'chn_sat', lf_dataparam.chn_sat, ...
          'nframes', lf_dataparam.nframes, ...
          'nchns',   lf_dataparam.nchns);
end

% Done
te_load = toc(ts_load);
fprintf('LUMO file loaded in %.1fs\n', te_load);

end



% load_lumo_chdat
%
% Load raw channel data intensity binary files.
function [chn_dat] = load_lumo_data(lf_dir, lf_desc, dataparams)

chn_dat = zeros(dataparams.nchns, dataparams.nframes, 'single');
chn_off = 1;

for fi = 1:length(lf_desc.intensity_desc)
  
  fid = fopen(fullfile(lf_dir, lf_desc.intensity_desc(fi).fn));
  
  magic = fread(fid, 1, '*uint8');
  binver = fread(fid, 3, 'uint8');
  fseek(fid, 4, 'cof');
  nchns = fread(fid, 1, 'uint64');
  nframes = fread(fid, 1, 'uint64');
  fseek(fid, 20, 'cof');
  binfin = ~fread(fid, 1, 'uint8=>logical');
  bigend = fread(fid, 1, 'uint8=>logical');
  fseek(fid, 2, 'cof');
  
  if magic ~= 0x92
    error('LUMO file (%s): intensity data contains invalid identifier in file %s\n', ...
      lf_dir, lf_desc.intensity_desc(fi).fn);
  end
  
  if ~(all(binver == [0x00 0x01 0x00].') || all(binver == [0x00 0x00 0x01].'))
    error('LUMO file (%s): intensity data of unknown version in file %s\n', ...
      lf_dir, lf_desc.intensity_desc(fi).fn);
  end
  
  
  if bigend
    error('LUMO file (%s): intensity data is big-endian in file %s\n', ...
      lf_dir, lf_desc.intensity_desc(fi).fn);
  end
  
  if (nchns ~= dataparams.nchns)
    error('LUMO file (%s): intensity data size mismatch in file %s\n', ...
      lf_dir, lf_desc.intensity_desc(fi).fn);
  end
  
  if (fi == length(lf_desc.intensity_desc))
    if(binfin ~= true)
      error('LUMO file (%s): intensity data final flag not set on last chunk in file %s\n', ...
        lf_dir, lf_desc.intensity_desc(fi).fn);
    end
  else
    if(binfin == true)
      error('LUMO file (%s): intensity data final flag set before last chunk in file %s\n', ...
        lf_dir, lf_desc.intensity_desc(fi).fn);
    end
  end
  
  if((nframes + chn_off - 1) > size(chn_dat,2))
    error('LUMO file (%s): intensity data exceeds reported frame count on file %s\n', ...
      lf_dir, lf_desc.intensity_desc(fi).fn);
  end
  
  chn_dat(:, chn_off:(chn_off+nframes-1)) = fread(fid, [nchns, nframes], '*single');
  chn_off = chn_off+nframes;
end

end


% lumo_load_events
%
% Construct events array from provided LUMO file.
function [events] = load_lumo_events(lf_dir, lf_desc)

% 1. Load and parse the event data
%
try
  raw = fileread(fullfile(lf_dir, lf_desc.ev_fn));
  events_raw = lumofile.toml.decode(raw);
catch e
  fprintf(2, 'LUMO file (%s) invalid: error parsing events file %s\n', lf_dir, lf_desc.ev_fn);
  rethrow(e);
end

% Check for an empty structure (no events)
if isempty(fieldnames(events_raw))
  events = [];
else
  
  lf_events = reqfield(events_raw, 'events', lf_dir);
  
  try
    for ii = 1:length(lf_events)
      
      if (lf_desc.lfver(1) == 0) && (lf_desc.lfver(2) < 4)
        % Versions prior to 0.4.0 stored the timestamp as a string (contrary to spec)
        ts = str2double(lf_events{ii}.Timestamp);
      else
        ts = lf_events{ii}.Timestamp;
      end

      % Note that timestamps are in seconds
      events(ii) = struct('mark', lf_events{ii}.name, 'timestamp', ts);
      
    end
  catch e
    fprintf(2, 'LUMO file (%s): error parsing events structure from file %s\n', lf_dir, lf_desc.ev_fn);
    rethrow(e);
  end
  
  % Sort by timestamp
  [~, perm] = sort([events.timestamp]);
  events = events(perm);
  
end

fprintf('LUMO file events file contains %d entries\n', length(events));

end





