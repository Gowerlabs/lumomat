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
%   'ignore_memory':  A logical value which indicates if the function should ignore any
%                     potential performance issues arising from loading large amounts of
%                     data. Defaults to false.
%
%   Returns:
%
%     enum:   An enumeration of the system containing:
%
%             enum.hub:     a description of the LUMO Hub used for recording
%             enum.groups:  an array of structures describing each group (cap) connected.
%                           LUMO files currently only store a single group, so this array
%                           should be of length 1. The form of this structure is described
%                           further below.
%
%     data:   An array of structures of data form each group in the enumeration.
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
%   nodes = enum.groups(gi).nodes;
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
%   layout = enum.groups(gi).layout;
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

ts_load = tic;

% Parse inputs
p = inputParser;
addParameter(p, 'ignore_memory', false, @islogical);
addParameter(p, 'layout', []);
parse(p, varargin{:});

ignore_memory = p.Results.ignore_memory;
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
    warning([...
      'The sepcified LUMO file does not contain an embedded layout file, and no layout '...
      'has been specified when calling this function. The returned enumeration will '...
      'lack layout information, and it will not be possible to convert this file to '...
      'formats which require a layout. Specify an appropriate layout file to supress '...
      'this warning, or copy an appropriate layout to the .LUMO folder in order for it '...
      'to be used as an automatic fallback.']);
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

% Load intensity files
%
% We first perform a check to ensure that the data structures to be created are not larger
% than the available system memory. This test can be overriden by the user, though evidently
% performance may suffer.
[~,sysmem] = memory;
mem_avail_mib = sysmem.PhysicalMemory.Available/1024/1024;
mem_reqrd_mib = lf_dataparam.nframes * lf_dataparam.nchns * 4 /1024/1024;

if  (mem_reqrd_mib < mem_avail_mib) || ignore_memory
  [chn_dat] = load_lumo_data(lf_dir, lf_desc, lf_dataparam);
else
  error([...
    'Loading intensity data for file %s will exceed available system memory. This error '...
    'be supressed by passing the argument pair (''ignore_memory'', true), but performance ' ...
    'may be impacted']);
end


data = struct('chn_dat', chn_dat, ...
  'chn_fps', lf_dataparam.chn_fps, ...
  'chn_dt',  round((1/lf_dataparam.chn_fps)*1000), ...
  'chn_sat', lf_dataparam.chn_sat, ...
  'nframes', lf_dataparam.nframes, ...
  'nchns',   lf_dataparam.nchns);

% Done
te_load = toc(ts_load);
fprintf('LUMO file loaded in %.1fs\n', te_load);

end


% load_lumo_desc
%
% Create the lumo file description, returning a structure of information required for
% parsing, including the version, file-names, and flags indicating the available contents.
% All files referenced in the structure have been confirmed to exist within the file.
%
function [lf_desc] = load_lumo_desc(lf_dir)

% 0. Some constants
lf_known_ver = [0 0 1; 0 1 0; 0 1 1; 0 2 0; 0 3 0; 0 4 0];

% 1. Ensure that the specified file exists
%
if exist(lf_dir, 'dir') ~= 7
  error('The specified LUMO file (%s) cannot be found', lf_dir)
else
  fprintf('Loading LUMO file %s\n', lf_dir);
end


% 2. Load and validate file metadata
%
% The metadata filenmae is fixed and if it cannot be found there is no way to proceed with
% automatic loading or conversation of the data.
%

% 2a. Attempt to load and parse the metadata
%
metadata_fn = fullfile(lf_dir, 'metadata.toml');
if exist(metadata_fn, 'file') ~= 2
  error('LUMO file (%s): invalid metadata not found\n', lf_dir)
end

try
  raw = fileread(metadata_fn);
  metadata = lumofile.toml.decode(raw);
catch e
  fprintf(2, 'LUMO file (%s): error parsing metadata file %s\n', lf_dir, metadata_fn);
  rethrow(e);
end

% 2b. Parse version information and check we understand this format
%
lf_ver = reqfield(metadata, 'lumo_file_version', lf_dir);
try
  lf_ver_num = str2double(strsplit(lf_ver, '.'));
catch e
  fprintf(2, 'LUMO file (%s): error parsing fiile version number\n', lf_dir);
  rethrow(e)
end

% We need to get the (required) file name field here, in order to deal with a version number
% ambiguity whereby 0.2.0 files can report as 0.1.0.
lf_meta_fns = reqfield(metadata, 'file_names');
if all(lf_ver_num == [0 1 1])
  if isfield(lf_meta_fns, 'hardware_file')
    lf_ver_num = [0 2 0];
  end
end

if ~ismember(lf_known_ver, lf_ver_num, 'rows')
  error('LUMO file (%s): version %d.%d.%d is not supported by this software version\n', ...
    lf_dir, lf_ver_num(1), lf_ver_num(2), lf_ver_num(3));
else
  fprintf('LUMO file version %d.%d.%d\n', ...
    lf_ver_num(1), lf_ver_num(2), lf_ver_num(3));
end

% 2c. Get filenames of additional metadata and check existence
%

% Look for the file which enuemrates the system
if (lf_ver_num(1) < 1) && (lf_ver_num(2) < 2)
  lf_meta_hw_fn = reqfield(lf_meta_fns, 'layout_file');
else
  lf_meta_hw_fn = reqfield(lf_meta_fns, 'hardware_file');
end

% Look for the file which describes the nature of the recording
if (lf_ver_num(1) < 1) && (lf_ver_num(2) < 1)
  % In versions <= 0.1.0 the 'recording data' file was known as the 'sd' file.
  lf_meta_hw_fn = reqfield(lf_meta_fns, 'sd_file');
else
  lf_meta_rd_fn = reqfield(lf_meta_fns, 'recordingdata_file');
end

% Look for log files
lf_meta_lg_fn = optfield(lf_meta_fns, 'log_file');
if isempty(lf_meta_lg_fn)
  fprintf('LUMO file %s does not contain log information\n', lf_dir);
  lf_has_log = false;
else
  lf_has_log = true;
end

% Look for a cap layout file
if (lf_ver_num(1) < 1) && (lf_ver_num(2) > 1)
  
  % On file versions >= 0.2.0 the layout file will be specified, or missing
  lf_meta_lo_fn = optfield(lf_meta_fns, 'layout_file');
  if isempty(lf_meta_lo_fn)
    fprintf('LUMO file does not contain layout information\n');
    lf_has_layout = false;
  else
    lf_has_layout = true;
  end
  
else
  
  % On file versions < 0.2.0 the layout file field is reserved for use by the enumeration
  % but we might still seek a layout by file name
  lf_meta_lo_fn = fullfile(lf_dir, 'layout.json');
  if exist(lf_meta_lo_fn, 'file') ~= 2
    fprintf('LUMO file %s does not contain cap layout information\n', lf_dir);
    lf_has_layout = false;
    lf_meta_lo_fn = [];
  else
    lf_has_layout = true;
  end
  
end

% Look for option event file
lf_meta_ev_fn = optfield(lf_meta_fns, 'event_file');
if isempty(lf_meta_ev_fn)
  fprintf('LUMO file %s does not contain event information\n', lf_dir);
  lf_has_events = false;
else
  lf_has_events = true;
end

% Make sure that all available files exist!
if exist(fullfile(lf_dir, lf_meta_rd_fn), 'file') ~= 2
  error('LUMO file (%s) invalid: recording file %s not found', lf_dir, lf_meta_rd_fn);
end

if exist(fullfile(lf_dir, lf_meta_hw_fn), 'file') ~= 2
  error('LUMO file (%s) invalid: hardware file %s not found', lf_dir, lf_meta_hw_fn);
end

if lf_has_log
  if exist(fullfile(lf_dir, lf_meta_lg_fn), 'file') ~= 2
    error('LUMO file (%s) invalid: log file %s not found', lf_dir, lf_meta_lg_fn);
  end
end

if lf_has_layout
  if exist(fullfile(lf_dir, lf_meta_lo_fn), 'file') ~= 2
    % Some layout files are saved with uppercase extension, which matters on some platforms
    [~, flname, flext] = fileparts(lf_meta_lo_fn);
    lf_meta_lo_fn = [flname upper(flext)];
    if exist(fullfile(lf_dir, lf_meta_lo_fn), 'file') ~= 2
      error('LUMO file (%s) invalid: layout file %s not found', lf_dir, lf_meta_lo_fn);
    end
  end
end

if lf_has_events
  if exist(fullfile(lf_dir, lf_meta_ev_fn), 'file') ~= 2
    error('LUMO file (%s) invalid: event file %s not found', lf_dir, lf_meta_ev_fn);
  end
end


% 3. Get filenames for intensity files, and check existence
lf_intensity_fns = optfield(metadata, 'intensity_files');

if isempty(lf_intensity_fns)
  error('LUMO file (%s) does not contain any channel intensity measurements', lf_dir);
end

try
  for ii = 1:length(lf_intensity_fns)
    
    lf_int_fn = lf_intensity_fns{ii}.file_name;
    lf_int_tr = lf_intensity_fns{ii}.time_range;
    
    if exist(fullfile(lf_dir, lf_int_fn), 'file') ~= 2
      error('LUMO file (%s) invalid: linked intensity file %s not found', lf_dir, lf_int_fn);
    end
    
    lf_int_ts(ii) = lf_int_tr(1);
    lf_intensity_desc(ii) = struct('fn', lf_int_fn, 'tr', lf_int_tr);
    
  end
catch e
  fprintf(2, 'LUMO file (%s) invalid: error parsing intensity structure from file %s\n', ...
    lf_dir, lf_meta_rd_fn);
  rethrow(e);
end

% Sort by start time
[~, perm] = sort(lf_int_ts);
lf_intensity_desc = lf_intensity_desc(perm);


% 4. Form output structure
lf_desc = struct('lfver', lf_ver_num, ...
  'rd_fn', lf_meta_rd_fn, ...
  'hw_fn', lf_meta_hw_fn, ...
  'lo_fn', lf_meta_lo_fn, ...
  'ev_fn', lf_meta_ev_fn, ...
  'lg_fn', lf_meta_lg_fn, ...
  'has_log', lf_has_log, ...
  'has_layout', lf_has_layout, ...
  'has_events', lf_has_events, ...
  'intensity_desc', lf_intensity_desc);

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




% load_lumo_enum
%
% Construct the canonical enumeration from the metadata in the specified LUMO file, in
% addition to parameters of the data which are contained in the same metadata files.
%
function [enum, dataparams] = load_lumo_enum(lf_dir, lf_desc)

% 1. Load and parse the enumeration data
%
try
  if (lf_desc.lfver(1) < 0) && (lf_desc.lfver(2) < 4)
    raw = fileread(fullfile(lf_dir, lf_desc.hw_fn));
    enum_raw = lumofile.toml.decode(raw);
  else
    enum_raw = toml_read_fixup_hw(fullfile(lf_dir, lf_desc.hw_fn));
  end
catch e
  fprintf(2, 'LUMO file (%s) invalid: error parsing (fixed) hardware file %s\n', lf_dir, lf_desc.hw_fn);
  rethrow(e);
end

% 2. Canonicalise the enumeration
%
enum = struct();

% 2a. Enumerate the hub
%
try
  
  if isfield(enum_raw.Hub, 'firmware_version')
    temp = enum_raw.Hub.firmware_version;
    if temp ~= ""
      enum.hub.fw_ver = string(temp);
    else
      enum.hub.fw_ver = '';
    end
  else
    enum.hub.fw_ver = '';
  end
  
  if isfield(enum_raw.Hub, 'hardware_version')
    temp = enum_raw.Hub.hardware_version;
    if temp ~= ""
      enum.hub.type = temp;
    else
      enum.hub.type = [];
    end
  else
    enum.hub.type = [];
  end
  
  if isfield(enum_raw.Hub, 'hub_serial_number')
    temp = enum_raw.Hub.hub_serial_number;
    if temp ~= -1
      enum.hub.serial = string(temp);
    else
      enum.hub.sn = '';
    end
  else
    enum.hub.sn = '';
  end
  
  
catch e
  fprintf(2, 'LUMO file (%s): an exception ocurred parsing the hub enumeration\n', lf_dir);
  rethrow(e);
end

% 2b. Enumerate the groups
%
% Since the LUMO file embeds the global indexing in the hardware
% description, we will extract the global indices in order to form a
% mapping for local <-> global.
%
try
  
  ng =  length(enum_raw.Hub.Group);
  if (ng < 1) || (ng > 1)
    error('LUMO file (%s) invalid: file contains more than one group', lf_dir);
  end
  
  % Over each group
  for gi = 1:length(enum_raw.Hub.Group)
    
    enum.groups(gi) = struct();
    
    % Normalise the UID to a hex string and a name
    [uid_hex, uid_name] = lumomat.norm_gid(enum_raw.Hub.Group(1).uid);
    enum.groups(gi).id = uid_hex;
    enum.groups(gi).name = uid_name;
    
    
    % Over each node (first pass for sorting, and accumulating source/detector count)
    nn = length(enum_raw.Hub.Group(gi).Node);
    nodeidsort = zeros(nn, 1);
    raw_gnq = 0;  % Total number of sources in raw enumeration
    raw_gnm = 0;  % Total number of detectors in raw enumeration
    raw_gnw = 2;  % TODO: wavelength number is fixed, this must be gaurded
    for ni = 1:nn
      nodeidsort(ni) = enum_raw.Hub.Group(gi).Node{ni}.node_id;
       raw_gnq = raw_gnq + length(enum_raw.Hub.Group(gi).Node{ni}.Source)/raw_gnw;
       raw_gnm = raw_gnm + length(enum_raw.Hub.Group(gi).Node{ni}.Detector);
    end
    [~, node_perm] = sort(nodeidsort);
    
    src_g2l = zeros(raw_gnq, raw_gnw, 2); % Global to local mapping for sources
    det_g2l = zeros(raw_gnm, 2);          % Global to local mapping for detectors
    src_pwr_g = zeros(raw_gnq, raw_gnw);
    
    % Over each node(second pass for construction)
    for nii = 1:nn
      
      ni = node_perm(nii);
      
      cs_node = struct(...
        'id',       enum_raw.Hub.Group(gi).Node{ni}.node_id, ...
        'rev',      enum_raw.Hub.Group(gi).Node{ni}.revision_id, ...
        'fwver',    enum_raw.Hub.Group(gi).Node{ni}.firmware_version);
      
      % Build array of sources
      lf_srcs = [enum_raw.Hub.Group(gi).Node{ni}.Source{:}];
      if length(lf_srcs) ~= 6
        error('LUMO file (%s): source enumeration error', lf_dir);
      end
      
      % Build the source list in the internal indexing layout, which differs from the
      % ordering used in the LUMO global-spectroscopic format.
      for qi = 1:length(lf_srcs)
        [cs_srcs_temp, src_idx, src_gi, src_gw, src_pwr] = trans_src_desc(lf_srcs(qi), lf_desc.lfver, lf_dir);
        cs_srcs(src_idx) = cs_srcs_temp;
                
        % Map from LUMO global-spectroscopic to node local
        %
        % This will subsequently be used to map the global system of channels into the node
        % local set.
        assert(all(src_g2l(src_gi, src_gw, :) == 0));
        src_g2l(src_gi, src_gw, :) = [nii src_idx];
        assert(src_pwr_g(src_gi, src_gw) == 0);
        src_pwr_g(src_gi, src_gw) = src_pwr;
      end
      
      % Build array of detectors
      lf_dets = [enum_raw.Hub.Group(gi).Node{ni}.Detector{:}];
      if length(lf_dets) ~= 4
        error('LUMO file (%s): detector enumeration error', lf_dir);
      end
      
      for mi = 1:length(lf_dets)
        [cs_dets_temp, det_idx, det_gi] = trans_det_desc(lf_dets(mi), lf_desc.lfver, lf_dir);
        cs_dets(det_idx) = cs_dets_temp;
        assert(all(det_g2l(det_gi, :) == 0));
        det_g2l(det_gi,:) = [nii det_idx];
      end
      
      % Build array of optodes.
      %
      % This infromation is hard coded, but checks are performed during the construction of
      % the source and detector arrays to ensure that this representation is consistent.
      % Additionally, we check here that the representation is complete.
      srcopt = [cs_srcs(:).optode_idx];
      detopt = [cs_dets(:).optode_idx];
      optidc = sort(unique([srcopt detopt]));
      if any(optidc ~= 1:7)
        error('LUMO file (%s): optode enumeration error', lf_dir);
      end
      for oi = 1:7
        cs_opts(oi) = gen_optode_desc(oi, lf_dir);
      end
      
      % Perform assignment
      cs_node.srcs = cs_srcs;
      cs_node.dets = cs_dets;
      cs_node.optodes = cs_opts;
      enum.groups(gi).nodes(ni) = cs_node;
      
    end
    
    % Assert sorted by ID
    [~, perm] = sort([enum.groups(gi).nodes.id]);
    assert(all(perm == 1:length(enum.groups(gi).nodes)));
    
  end
  
catch e
  fprintf(2, 'LUMO file (%s): an exception ocurred parsing the group enumeration\n', lf_dir);
  rethrow(e);
end


% 2c. Enumerate the channels
%
% The LUMO file enumerates the channel data using global indices, and we seek the canonical
% node local representation. We must thus map back from the global -> local mapping.
%
try
  if (lf_desc.lfver(1) < 1) && (lf_desc.lfver(2) < 4)
    raw = fileread(fullfile(lf_dir, lf_desc.rd_fn));
    rcdata = lumofile.toml.decode(raw);
  else
    rcdata = toml_read_fixup_rc(fullfile(lf_dir, lf_desc.rd_fn));
  end
catch e
  fprintf(2, 'LUMO file (%s): error parsing recording data file %s\n', lf_dir, lf_desc.rd_fn);
  rethrow(e);
end

% Note that from here on in we are certainly using group index 0 (1) as there is presently
% no way to encode more than one group in a .lumo
gi = 1;

% Assert redundant data
try
  
  lf_nodes = rcdata.variables.nodes;
  lf_nsrc = rcdata.variables.n_srcs;
  lf_ndet = rcdata.variables.n_dets;
  lf_wls = rcdata.variables.wavelength;
  
  assert(all(sort(lf_nodes) == [enum.groups(gi).nodes.id]));
  assert(lf_nsrc == size(src_g2l, 1));
  assert(lf_ndet == size(det_g2l, 1));
  
  for ni = 1:length(enum.groups(gi).nodes)
    assert(unique(all(lf_wls == sort(unique([enum.groups(gi).nodes(1).srcs.wl])))));
  end
  
  % This assert is required to ensure the global to local indexing is valid
  assert(all(lf_wls == [735 850]))
  
catch e
  fprintf(2, 'LUMO file (%s): an exception ocurred whilst asserting group consistency\n', lf_dir);
  rethrow(e);
end

% Process the channel list and form the canonical enumeration
try
  
  lf_nchans = rcdata.variables.n_chans;
  lf_chlist = rcdata.variables.chans_list; %reshape(, 3, lf_nchans).';
  lf_chvalid = rcdata.variables.chans_list_act.';
  lf_t0 = rcdata.variables.t_0;
  lf_tlast = rcdata.variables.t_last;
  lf_framerate = rcdata.variables.framerate;
  lf_nframes = rcdata.variables.number_of_frames;
  
  assert(size(lf_chvalid, 1) == lf_nchans);
  assert(size(lf_chlist, 1) == lf_nchans);
  assert(size(lf_chvalid, 1) == lf_nchans);
  
catch e
  fprintf(2, 'LUMO file (%s): an exception ocurred whilst building channel enumeration\n', lf_dir);
  rethrow(e);
end

for ci = 1:lf_nchans
  
  % Get the global indices
  g_srcidx = lf_chlist(ci, 1);
  g_detidx = lf_chlist(ci, 2);
  g_wlidx = lf_chlist(ci, 3);
  
  % Map global to local sources
  l_src_node_idx = src_g2l(g_srcidx, g_wlidx, 1);
  l_src_node_id = enum.groups(gi).nodes(l_src_node_idx).id;
  l_src_idx = src_g2l(g_srcidx, g_wlidx, 2);
  l_src_optode_idx = enum.groups(gi).nodes(l_src_node_idx).srcs(l_src_idx).optode_idx;
  
  % Map global to local detectors
  l_det_node_idx = det_g2l(g_detidx, 1);
  l_det_node_id = enum.groups(gi).nodes(l_det_node_idx).id;
  l_det_idx = det_g2l(g_detidx,2);
  l_det_optode_idx = enum.groups(gi).nodes(l_det_node_idx).dets(l_det_idx).optode_idx;
  
  % Assign the node local and global-spectroscopic indices
  cs_chan = struct(...
    'src_node_idx', l_src_node_idx, ...
    'src_idx', l_src_idx, ...
    'det_node_idx', l_det_node_idx, ...
    'det_idx', l_det_idx);
  
  
  % The following fields are useful when working data, but they are derived paramters and
  % can be found by indexing the enumeration quite easily.
  %
  % 'src_node_id', l_src_node_id, ...
  % 'src_optode_idx', l_src_optode_idx,...
  % 'det_node_id', l_det_node_id, ...
  % 'det_optode_idx', l_det_optode_idx
  
  enum.groups(gi).channels(ci) = cs_chan;
  
end

if isfield(enum.groups(gi),'layout')
  lumomat.validate_layout(enum);
end

fprintf('LUMO file enumeration contains %d tiles, %d channels\n', ...
  length(enum.groups(gi).nodes), length(enum.groups(gi).channels));


dataparams = struct('nchns', lf_nchans, ...
  'nframes', lf_nframes, ...
  'chn_list', lf_chlist, ...
  'chn_sat', lf_chvalid, ...
  'chn_fps', lf_framerate);

end

% gen_optode_desc
%
% Generate a hard coded optode list for the optode index specified. This is consistent with
% the internal ordering of the LUMO stack.
function cs_opt = gen_optode_desc(io_idx, lf_dir)

switch io_idx
  case 1
    optode_name = '1';
    optode_type = 'D';
  case 2
    optode_name = '2';
    optode_type = 'D';
  case 3
    optode_name = '3';
    optode_type = 'D';
  case 4
    optode_name = '4';
    optode_type = 'D';
  case 5
    optode_name = 'A';
    optode_type = 'S';
  case 6
    optode_name = 'B';
    optode_type = 'S';
  case 7
    optode_name = 'C';
    optode_type = 'S';
  otherwise
    error('LUMO file (%s): error generating optode structure (opt id %d)',...
      lf_dir, io_idx);
end

cs_opt = struct('name', optode_name, 'type', optode_type);

end


% trans_src_desc
%
% Translate the .lumo file source descriptions into the canonical node local format.
% Assertions are required to allow us to build a hard-coded optode array.
%
% Function outputs the internal local source index, the reported global-spectroscopic source
% index, wavelength index, and source power.
%
% Notes:
%
%   - The wavelength indexing (src_gwi) is hard coded, it must be asserted
%     that the ordering (1 -> 735, 2 -> 850) is valid.
%
function [cs_src, src_idx, src_gsi, src_gwi, src_pwr] = trans_src_desc(lf_src, lf_ver, lf_dir)

% ID depends upon metadata version
if (lf_ver(1) < 1) && (lf_ver(2) < 4)
  idset = [1 2 4 8 16 32];
else
  idset = [1 1 5 5 3 3];
end

switch lf_src.description
  case 'SRCA_735nm'
    src_idx = 1;
    src_wl = 735;
    src_gwi = 1;
    src_optode_idx = 5;
    assert(lf_src.id == idset(1));
  case 'SRCA_850nm'
    src_idx = 4;
    src_wl = 850;
    src_gwi = 2;
    src_optode_idx = 5;
    assert(lf_src.id == idset(2));
  case 'SRCB_735nm'
    src_idx = 2;
    src_wl = 735;
    src_gwi = 1;
    src_optode_idx = 6;
    assert(lf_src.id == idset(3));
  case 'SRCB_850nm'
    src_idx = 5;
    src_wl = 850;
    src_gwi = 2;
    src_optode_idx = 6;
    assert(lf_src.id == idset(4));
  case 'SRCC_735nm'
    src_idx = 3;
    src_wl = 735;
    src_gwi = 1;
    src_optode_idx = 7;
    assert(lf_src.id == idset(5));
  case 'SRCC_850nm'
    src_idx = 6;
    src_wl = 850;
    src_gwi = 2;
    src_optode_idx = 7;
    assert(lf_src.id == idset(6));
  otherwise
    error('LUMO file (%s): error parsing source structure (src id %d)', lf_dir, lf_src.id);
end

src_pwr = lf_src.Source_power;
src_gsi = lf_src.group_location_index;
cs_src = struct('wl', src_wl, 'optode_idx', src_optode_idx, 'power', src_pwr);

end


% trans_det_desc
%
% Translate the .lumo detector descriptions into the canonical node local format. All
% indices are zero based. Assertions are required to allow us to build a hard-coded optode
% array.
%
% Also outputs the global detector index.
function [cs_det, det_idx, det_gsi] = trans_det_desc(lf_det, lf_ver, lf_dir)


adcbase = 'ADC Detector Channel ';

% ID depends upon metadata version
if (lf_ver(1) < 1) && (lf_ver(2) < 4)
  adcdesc = {[adcbase '0'], [adcbase '1'], [adcbase '2'], [adcbase '3']};
  idset = [0 1 2 3];
else
  adcdesc = {[adcbase '1'], [adcbase '2'], [adcbase '3'], [adcbase '4']};
  idset = [4 6 0 2];
end

switch lf_det.description
  case adcdesc(1)
    det_idx = 1;
    det_optode_idx = 1;
    assert(lf_det.id == idset(1));
  case adcdesc(2)
    det_idx = 2;
    det_optode_idx = 2;
    assert(lf_det.id == idset(2));
  case adcdesc(3)
    det_idx = 3;
    det_optode_idx = 3;
    assert(lf_det.id == idset(3));
  case adcdesc(4)
    det_idx = 4;
    det_optode_idx = 4;
    assert(lf_det.id == idset(4));
  otherwise
    error('LUMO file (%s): error parsing detector structure (det id %d)', lf_dir, lf_det.id);
end

det_gsi = lf_det.group_location_index;
cs_det = struct('optode_idx', det_optode_idx);

end


% reqfield
%
% Get a required field from a structure, throwing an appropriate error if unavailable.
%
function fieldval = reqfield(s, name, lumodir_fn)

if isfield(s, name)
  fieldval = getfield(s, name);
else
  error('LUMO file (%s): required field %s missing from structure', lumodir_fn, name)
end

end


% reqfieldci
%
% Get a required field from a structure, case insensitive on the field name, ithrowing an
% appropriate error if unavailable.
%
function fieldval = reqfieldci(s, name, lumodir_fn)

names   = fieldnames(s);
isField = strcmpi(name,names);

if any(isField)
  fieldval = s.(names{isField});
else
  fieldval = [];
end

if isempty(fieldval)
  error('LUMO file (%s): required field %s missing from structure', lumodir_fn, name)
end

end

% optfield
%
% Attempt to get an optional field from a structure, returning an empty array if unavailable.
%
function fieldval = optfield(s, name)

if isfield(s, name)
  fieldval = getfield(s, name);
else
  fieldval = [];
end

end


% toml_fixup_hardware
%
% Modify whitespace in toml file to ensure compatibility with matlab-toml parser.
%
function toml_data = toml_read_fixup_hw(fn)

raw_text = fileread(fn);
fixed_raw_text = regexprep(raw_text,'[\n\r]+[\t ]+','\n');
toml_data = lumofile.toml.decode(fixed_raw_text);

end

% toml_read_fixup_rc
%
% Modify whiltespace in toml file to ensure compatiblity with matlab-toml parse. This
% function is typically applied to the recording data files.
%
function toml_data = toml_read_fixup_rc(fn)

raw_text = fileread(fn);

% Fix arrays with new lines and indentation before numbers or opening brackets.
raw_text = regexprep(raw_text,'[\n\r]+[\t ]+([\d\[])','$1');

% Fix arrays with newlines before closing bracket.
raw_text = regexprep(raw_text,'[\n\r]+\]',']');

% Removes spaces between delimiters and numeric elements of an array.
raw_text = regexprep(raw_text,'([^=]) +(\d+)','$1$2');

% Removes spaces between numeric elements of an array and delimters.
raw_text = regexprep(raw_text,'(\d+) +([^=])','$1$2');

% Adds spaces between, delimiters and the numeric element that comes after.
raw_text = regexprep(raw_text,',(.)',', $1');

toml_data = lumofile.toml.decode(raw_text);

end



