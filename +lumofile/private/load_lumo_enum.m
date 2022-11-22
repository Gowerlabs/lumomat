function [enum, dataparams] = load_lumo_enum(lf_dir, lf_desc)
% LOAD_LUMO_ENUM
%
% Construct the canonical enumeration from the metadata in the specified LUMO file, in
% addition to parameters of the data which are contained in the same metadata files.
%

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

% LUMO file versions 0.4.0&0.5.0 don't have channel saturation data
if all(lf_desc.lfver == [0 4 0]) || all(lf_desc.lfver == [0 5 0])
    fprintf('LUMO file does not contain saturation data\n');
    dataparams = struct('nchns', lf_nchans, ...
      'nframes', lf_nframes, ...
      'chn_list', lf_chlist, ...
      'chn_fps', lf_framerate);
else
    dataparams = struct('nchns', lf_nchans, ...
      'nframes', lf_nframes, ...
      'chn_list', lf_chlist, ...
      'chn_sat', ~lf_chvalid, ...
      'chn_fps', lf_framerate);
end

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