function [enum, data, events] = read_lufr(lufrfn, varargin)
% LUMOFILE.READ_LUFR Read a LUFR file from disk
%
% [enum, data, events] = LUMOFILE.READ_LUMO(filename) 
%
% LUMOFILE.READ_LUFR is used to read data the LUFR file format emitted by Gowerlabs
% development tools for LUMO. 
%
%   Parameters
%
%   fn: input filename
%
%   Optional parameters
%
%   'group':     LUFR files may store more than one group ('cap'), but this script will
%                only output the raw data from a single group at a time. Choose the index 
%                using this  variable.
%
%
%   'layout':    When a layout is specified, it is embedded into the output enumeration for
%                direct use or for export. A layout can be specified as either:
%
%                       string: the filename of a valid lumo layout file, in JSON format
%                       struct: a layout structure in the format returned by
%                               lumofile.read_layout (see function help for details)
%
%   'optfilter': When loading the data, keep only those channels which are recorded on 
%                optode pairs provided in the specified matrix. The matrix should have rows 
%                of the form:
%
%                [ src_node_id src_opt_id det_node_id det_opt_id]
%                   
%                where:
%                 - src_node_id is the one-indexed dock ID of the source node
%                 - src_opt_id is the source optode, indexed as A=1, B=2, C=3
%                 - det_node_id is the one-indexed dock ID of the detector node
%                 - det_opt_id is the detector optode index, [1, 4].
%
%   Returns:
%
%     enum:   An enumeration of the system containing:
%
%             enum.hub:     a description of the LUMO Hub used for recording
%             enum.groups:  an array of structures describing each group (cap) connected. Since
%                           we will only load a single group at a time, this structure is scalar
%                           and indexing is not reuiqred.
%
%     data:   A structure of data form the selected group.
%
%     events: An array of structures detailing events recorded during recording.
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
%   >> ch = enum.groups.channels(98)
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
%   The data:
%
%   The primary purpose of the data structure is to provide the channel intensity data and
%   fixed associated metadata.
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

% Normalise strings
lufrfn = convertStringsToChars(lufrfn);

% Parse the inputs
apply_filter = false;

p = inputParser;
addParameter(p, 'layout', []);
addParameter(p, 'group', 1, @isnumeric);
addParameter(p, 'chfilter', [], @ismatrix);
parse(p, varargin{:});

chfilter = p.Results.chfilter;       % Channel optode-pair filtering matrix
gidx = p.Results.group-1;            % Selected group index for hyperscanning data
layout_override = p.Results.layout;

if exist(lufrfn, 'file') ~= 2
  error('The specified LUFR file (%s) cannot be found', lufrfn)
else
  fprintf('Loading LUFR file %s\n', lufrfn);
end


if(~isempty(chfilter))
    
    if(size(chfilter,2) ~= 4)
        error('Proivded channel filter matrix is invalid');
    end
    
    apply_filter = true;
    
    fprintf('Applying optode filter with %d entries\n', size(chfilter,1));
    
end

% Assorted constants
tag_information = 1;
tag_enumeration = 2;
tag_frame = 3;
tag_event = 4;
tag_logentry = 5;   % Introduced in v3
tag_layout = 6;     % Introduced in v3
tag_max = tag_layout;

% Record the minimum number of wavelengths expected, which is used
% during channel filtering.
n_wl_min = 2;

% Open the file
fid = fopen(lufrfn, 'rb');

% Get file size
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

% Check the file is valid
filehdr = string(fread(fid, 2, 'int8=>char').');
if(filehdr ~= "LU")
    error('The lumo frame file does not contain the correct file header');
end

% Check the file version
filever = fread(fid, 1, 'uint8=>double');
if(filever ~= 1 && filever ~=2 && filever ~=3)
    error('The lumo frame file is of an unknown version %i', filever);
else
  fprintf('LUFR file version %d\n', filever);
end

% On version two, we will need to set the particular group index from
% which we wish to extract data
if(filever > 1)
    fprintf('LUFR file group index %d selected \n', gidx);
end

% Skip over zeros
if(filever > 2)
  flags = fread(fid, 3, 'uint8=>double');
  if(flags(2) == 1)
    warning('Gain lock low enabled');
  end
  if(flags(2) == 2)
    warning('Gain lock high enabled');
  end
  
else
  fseek(fid, 3, 'cof');
end

% Get endieness marker
endmarker = fread(fid, 1, 'uint16=>uint16');
if(endmarker ~= hex2dec('1234'))
    error('The lumo frame file endienness is not supported');
end

% Read the record counter
record_count = fread(fid, 1, 'uint32=>double');
if(record_count == 0)
    warning('Record counter not set, file may be corrupted');
end

% Now, seek through the file to count the number of frames, storing
% the offset in an array alongside the length
rcoffset = [];
rclength = [];
rctag = [];

% Also enumerate the number of particular record types
n_frames = 0;
n_enums = 0;
n_events = 0;
n_infos = 0;
n_logentries = 0;
n_layouts = 0;

% Special case the enumeration so we can access it easily
rc_enum_idx = 0;

i = 0;

while(~feof(fid))
    
    % Check for a tag
    recordtag = fread(fid, 1, 'uint32=>double');
    if(feof(fid))
        break;
    end
    
    if(recordtag > tag_max)
        if(feof(fid))
            break;
        end
        error('The lumo frame file contains a bad frame header');
    end
    
    % Get the length
    recordlen = fread(fid, 1, 'uint32=>double');
    
    % Record it
    rcoffset(end+1) = ftell(fid);
    rclength(end+1) = recordlen;
    rctag(end+1) = recordtag;
    
    % Special case the enumeration because we may need it early
    % for filtering. Note that we will reject files with multiple
    % enumerations so it doesn't matter if this is overwritten
    
    if(recordtag == tag_enumeration)
        fprintf('LUFR file enumeration found at tag %d\n', rc_enum_idx);
        rc_enum_idx = length(rctag);
    end
    
    % Check that standard frame record lengths are constant
    if(recordtag == tag_frame)
        
        % On file versions beyond 1, a group index is inserted into
        % frame records. If the group index is not that which has been
        % selected, we skip over the data
        if(filever > 1)
            group_idx = fread(fid, 1, 'int32=>doble');
            if(group_idx ~= gidx)
                fseek(fid, rcoffset(end) + recordlen, 'bof');
                continue;
            end
        end
        
        % On file versions beyond 1, an error counter is inserted between the frame
        % index and the remaining metadata.
        if(filever > 1)
            sizeparam = fread(fid, 10, 'int32=>double');
            sizeparam = sizeparam([1 3 4 5 6 7 8 9 10]);
        else
            sizeparam = fread(fid, 9, 'int32=>double');
        end
        
        if n_frames == 0
            % On the first run, we store this as a reference
            sizeparam_ref = sizeparam;
            firstframe = false;
        else
            % Subsequently, check no changes have happned
            if ~all(sizeparam(2:end) == sizeparam_ref(2:end))
                error('LUFR file invalid: frame data dimensions vary');
            end
            
        end
        
        n_frames = n_frames + 1;
        
    end
    
    if(recordtag == tag_enumeration)
        n_enums = n_enums + 1;
    end
    
    if(recordtag == tag_information)
        n_infos = n_infos + 1;
    end
    
    if(recordtag == tag_event)
        n_events = n_events + 1;
    end
    
    if(recordtag == tag_logentry)
      if(filever < 3)
          error('LUFR file of version < 3 contains a log entry record');
      end  
      n_logentries = n_logentries + 1;
    end
    
    if(recordtag == tag_layout)
        if(filever < 3)
          error('LUFR file of version < 3 contains a layout record');
        end
        
        group_idx = fread(fid, 1, 'int32=>doble');
          if(group_idx ~= gidx)
              fseek(fid, rcoffset(end) + recordlen, 'bof');
              continue;
          end
        n_layouts = n_layouts + 1;
    end
        
    % Skip to the next
    fseek(fid, rcoffset(end) + recordlen, 'bof');
    
    % Increment counter
    i = i+1;
 
end

if(n_enums < 1)
    error('LUFR file does not contain an enumeration block');
end

if(n_enums > 1)
    error('LUFR file contains more than one enumeration block');
end

if(n_layouts > 1)
    error('LUFR file contains more than one layout file for group index %i', group_idx);
end

if length(rclength) == record_count
    fprintf('LUFR file contains %d records, %d frames, %d events, %d layouts \n', length(rclength), n_frames, n_events, n_layouts);
else
    warning('LUFR file found %d records, %d frames, %d events in the input data file\n', length(rclength), n_frames, n_events);
end

% Work with the metadata

% Allocate information blocks
infoblks = strings(n_infos, 1);

% Allocate event times and strings
evtim = zeros(n_events, 1);
evstr = cell(n_events, 1);

% Get the enumeration
fseek(fid, rcoffset(rc_enum_idx), 'bof');

fprintf('Loading enumeration... ');
if(rctag(rc_enum_idx) == tag_enumeration)
    enumjson = fread(fid, rclength(rc_enum_idx), 'char=>char');
    enum = jsondecode(convertCharsToStrings(enumjson));
else
    fprintf('\n');
    error('Logical error loading enumeration data');
end

fprintf('done\n');


% Fix enmeration acquistion data in some versions
if(isfield(enum, 'group'))
    
    fprintf('LUFR file contains imperfect structure, fixing...\n');
    
    % Merge fields into correct structure. Note that versions with this
    % erroneous field naming did not have hyperscanning support, so the
    % group index can be fixed.
    for i = 1:length(enum.groups(1).channels)
        enum.groups(1).channels(i).acq_row = enum.group.channels(i).acq_row;
        enum.groups(1).channels(i).acq_offset_us = enum.group.channels(i).acq_offset_us;
    end
    
    enum = rmfield(enum, 'group');
    
end


% If we are to apply a channel filter, build the permutation maitrx now
if(apply_filter)
    
    n_filter = size(chfilter, 1);
    chperm = zeros(n_filter*n_wl_min, 1);
    n_schans_keep = 0;
    
    % Build the indexing arrays
    lin_src_node_id = [enum.groups(gidx + 1).channels.src_node_id].';
    lin_det_node_id = [enum.groups(gidx + 1).channels.det_node_id].';
    
    lin_src_opt = [enum.groups(gidx + 1).channels.src_optode_name].' - 64;    % ASCII 'A' -> 1
    
    if(filever > 2)
      lin_det_opt = [enum.groups(gidx + 1).channels.det_optode_name].' - 48;   % ASCII '1' -> 1
    else
      lin_det_opt = [enum.groups(gidx + 1).channels.det_optode_name].' - 47;   % ASCII '0' -> 1
    end
     
    % Over every entry in the keep list
    k = 1;
    for i = 1:size(chfilter,1)
                
        src_node_match = lin_src_node_id == chfilter(i,1);
        det_node_match = lin_det_node_id == chfilter(i,3);
        src_opt_match = lin_src_opt == chfilter(i,2);
        det_opt_match = lin_det_opt == chfilter(i,4);
        
        ch_match = find(src_node_match & det_node_match & src_opt_match & det_opt_match);
        n_ch_match = length(ch_match);
        
        chperm(k:(k+n_ch_match-1)) = ch_match;
        k = k+n_ch_match;
       
    end
    

    
    if(any(chperm == 0))
        error('Some entries in the channel keep filter could not be matched');
    end
    
    n_schans_keep = length(chperm);
    fprintf('LUFR file channel filtering complete, found %d channels\n', n_schans_keep);
    
else
    
    chperm = 1:length(enum.groups(gidx + 1).channels);
    
end

% Compute size of the output data an allocate
n_nodes  = sizeparam_ref(3);
n_schans = sizeparam_ref(4);
n_dchans = sizeparam_ref(5);
n_mpu    = sizeparam_ref(6);
n_spw    = sizeparam_ref(7);
n_det    = sizeparam_ref(8);
n_row    = sizeparam_ref(9);

% Get frames per second
tchdat = sizeparam_ref(2)*1e-6; %in seconds
fps = 1/(sizeparam_ref(2)*1e-6);

fprintf('LUFR file channel frame rate is %.2f fps\n', fps);

if(apply_filter)
    chdat = zeros(n_schans_keep, n_frames, 'single');   % Channel data
else
    chdat = zeros(n_schans, n_frames, 'single');        % Channel data
end
dkdat = zeros(n_dchans, n_frames, 'single');            % Dark channel data
tmpdat = zeros(n_nodes, n_frames, 'single');            % Temperature
vindat = zeros(n_nodes, n_frames, 'single');            % Input voltage
srcpwr = zeros(n_nodes, n_spw, n_frames, 'single');     % Source powers
detmax = zeros(n_nodes, n_row, n_det, n_frames, 'single'); % Detector maximum
accdat = zeros(n_nodes, 3, n_mpu, n_frames);            % MPU data (concatenated later)
gyrdat = zeros(n_nodes, 3, n_mpu, n_frames);            % MPU data (concatenated later)

errcnt = zeros(n_frames,1);                             % Errors in frame counter

% % Loop over the frames and get them data

i_fr = 1;   % Frame count
i_if = 1;   % Information block count
i_ev = 1;   % Event count

for i = 1:length(rclength)
    
    % Seek to the frame
    fseek(fid, rcoffset(i), 'bof');
    
    % Standard measurement frame
    if rctag(i) == tag_frame
        
        
        % Skip the record if it is not from the desired group
        if(filever > 1)
            group_idx = fread(fid, 1, 'int32=>doble');
            if(group_idx ~= gidx)
                continue;
            end
            
            % Get frame dimension data to locate the error counter
            sizeparam = fread(fid, 10, 'int32=>double');
            errcnt(i_fr) = sizeparam(2);
        else
            
            % Skip over the frame dimension data
            fread(fid, 9, 'int32=>double');
            
        end
        
        % Get the data
        if(apply_filter)
            chdat_temp = fread(fid, n_schans, 'single=>single');
            chdat(:,i_fr) = chdat_temp(chperm);
        else
            chdat(:,i_fr) = fread(fid, n_schans, 'single=>single');
        end
        
        dkdat(:,i_fr) = fread(fid, n_dchans, 'single=>single');
        
        % Skip over the MPU for now
        mpuframe = fread(fid, n_mpu*9*n_nodes, 'single=>single');
        mpudat_raw = permute(reshape(mpuframe, 9, n_mpu, n_nodes), [3 1 2]);
        
        for j = 1:n_mpu
            % mpudat_raw : [n_nodes, 9, n_mpu]
            gyrdat(:, :, j, i_fr) = mpudat_raw(:, 1:3, j);
            accdat(:, :, j, i_fr) = mpudat_raw(:, 4:6, j);
        end
        
        % Gerrabit more
        tmpdat(:,i_fr) = fread(fid, n_nodes, 'single=>single');
        vindat(:,i_fr) = fread(fid, n_nodes, 'single=>single');
        
        % Source powers
        srcpwr(:,:,i_fr) = reshape(fread(fid, n_spw*n_nodes, 'single=>single'), n_spw, n_nodes).';
        
        % DX maximum per detector on each frame
        detmax_i = reshape(fread(fid, n_det * n_row * n_nodes, 'single=>single'), n_det, n_row, n_nodes);
        detmax(:,:,:,i_fr) = permute(detmax_i, [3,2,1]);
        
        
        i_fr = i_fr+1;
    end
    
    if(rctag(i) == tag_event)
        evtim(i_ev) = fread(fid, 1, 'uint32=>double');
        evtxt = fread(fid, rclength(i)-4, 'char=>char');
        evstr{i_ev} = convertCharsToStrings(evtxt);
        i_ev = i_ev + 1;
        
    end
    
    if(rctag(i) == tag_information)
        infotxt = fread(fid, rclength(i), 'char=>char');
        infoblks(i_if) = convertCharsToStrings(infotxt);
        i_if = i_if + 1;
    end
    
    if(rctag(i) == tag_layout)
     
      group_idx = fread(fid, 1, 'int32=>doble');
      if(group_idx ~= gidx)
        continue;
      end
      layoutstr = fread(fid, rclength(i)-4, 'char=>char');
      
    end
    
end


% Close the file
fclose(fid);


% Build the saturation flag matrix
det_sat_limit = 97.5;

fprintf('LUFR file building saturation channel mapping... ');

% Build a [src_node_id x src x det_node_id x row x det] -> channel mapping
n_srcdim_max = 6;
satmap = zeros(size(detmax,1), n_srcdim_max, size(detmax,1), size(detmax,2), size(detmax,3));

for i = 1:length(enum.groups(gidx + 1).channels)
    i_row_idx = enum.groups(gidx+1).channels(i).acq_row + 1;
    i_det_node_idx = enum.groups(gidx+1).channels(i).det_node_idx + 1;
    i_det_opt_idx = enum.groups(gidx+1).channels(i).det_optode_idx + 1;
    i_src_node_idx = enum.groups(gidx+1).channels(i).src_node_idx + 1;
    i_src_opt_idx = enum.groups(gidx+1).channels(i).src_optode_idx + 1;
    if(enum.groups(gidx+1).channels(i).src_wl == 850)
        i_src_opt_idx = i_src_opt_idx-3;
    end
                    
    satmap(i_src_node_idx, i_src_opt_idx, i_det_node_idx, i_row_idx, i_det_opt_idx) = i;
end

if(apply_filter)    
     satflag = false(n_schans_keep, n_frames);
else
     satflag = false(n_schans, n_frames);
end   
                  
for j = 1:n_frames
    satindex = detmax(:,:,:,j) > det_sat_limit;
    satchans = satmap(:,:,satindex);
    satchans = satchans(satchans > 0);
    
    if(apply_filter)
        satflag_temp = false(n_schans, 1);
        satflag_temp(satchans) = true;
        satflag(:,j) = satflag_temp(chperm);
    else
        satflag(satchans,j) = true;
    end
    
end


fprintf('done\n');


% Reorder MPU data
dim = size(gyrdat);
gyrdat = reshape(gyrdat,[dim(1) dim(2) dim(3)*dim(4)]);
dim = size(accdat);
accdat = reshape(accdat,[dim(1) dim(2) dim(3)*dim(4)]);
tmpudat = 10e-3; %Always at 100Hz = 1/10ms


%%% Restructure the enumeration
%
% The loaded enumeration is direct from the stack and must be reorganised and clenaed to 
% conform to the canonical format.
%

% Rework the group ID to a hex string and descriptive name
[uid_hex, uid_name] = lumomat.norm_gid(enum.groups(gidx + 1).id);
enum.groups(gidx + 1).id = uid_hex;
if(filever > 2)
    fprintf("Embedded group name for index %u: %s\n", gidx, enum.groups(gidx + 1).name);
else
enum.groups(gidx + 1).name = uid_name;
end


% Convert all indices to one-based for MATLAB
nc = length(enum.groups(gidx + 1).channels);
for ci = 1:nc
  enum.groups(gidx + 1).channels(ci).det_idx = enum.groups(gidx + 1).channels(ci).det_idx + 1;
  enum.groups(gidx + 1).channels(ci).det_node_idx = enum.groups(gidx + 1).channels(ci).det_node_idx + 1;
  enum.groups(gidx + 1).channels(ci).src_idx = enum.groups(gidx + 1).channels(ci).src_idx + 1;
  enum.groups(gidx + 1).channels(ci).src_node_idx = enum.groups(gidx + 1).channels(ci).src_node_idx + 1;
  
  % On versions <= 2 there is a bug in which the optode indices are incorrect in the node
  % table, but they are correct in the channel table, we can replace one with 'tother and also
  % make one-based at the same time.
  if(filever < 3)
    enum.groups(gidx + 1).nodes(enum.groups(gidx + 1).channels(ci).src_node_idx).srcs(enum.groups(gidx + 1).channels(ci).src_idx).optode_idx =  enum.groups(gidx + 1).channels(ci).src_optode_idx;
    enum.groups(gidx + 1).nodes(enum.groups(gidx + 1).channels(ci).det_node_idx).dets(enum.groups(gidx + 1).channels(ci).det_idx).optode_idx =  enum.groups(gidx + 1).channels(ci).det_optode_idx;
  end
    
end

% Remove extraneous non-canonical fields from the channels structure
enum.groups(gidx + 1).channels = rmfield(enum.groups(gidx + 1).channels, ...
  {'acq_offset_us', 'acq_row',...
   'det_node_id', 'det_optode_idx', 'det_optode_name', ...
   'src_node_id', 'src_optode_idx', 'src_optode_name', 'src_wl'});

% Modify the nodes
nd = length(enum.groups(gidx + 1).nodes);
for ni = 1:nd
  % Convert node firmware revision to semver string, rename rev -> revision,
  fwmajor = enum.groups(gidx +1).nodes(ni).fwmajor;
  fwminor = enum.groups(gidx +1).nodes(ni).fwminor;
  fwpatch = enum.groups(gidx +1).nodes(ni).fwpatch;
  enum.groups(gidx + 1).nodes(ni).fwver = sprintf('%d.%d.%d', fwmajor, fwminor, fwpatch);
  
  % Remove extraneous fields
  if(filever < 3)
    enum.groups(gidx + 1).nodes(ni).optodes = rmfield(enum.groups(gidx + 1).nodes(ni).optodes, {'rho', 'theta'});
  end
  
  % Change optode naming from stack '0', '1', '2', '3', to canonical '1', '2', '3', '4', on
  % file versions < 3. On version >= 3 this is automatically done in the front end.
  no = length(enum.groups(gidx + 1).nodes(ni).optodes);
  for j = 1:no
    if(filever < 3)
      switch enum.groups(gidx + 1).nodes(ni).optodes(j).name
        case '0'
          enum.groups(gidx + 1).nodes(ni).optodes(j).name = '1';
        case '1'
          enum.groups(gidx + 1).nodes(ni).optodes(j).name = '2';
        case '2'
          enum.groups(gidx + 1).nodes(ni).optodes(j).name = '3';
        case '3'
          enum.groups(gidx + 1).nodes(ni).optodes(j).name = '4';
      end
    end    
  end

  % One base the the src and detector optode indices
  nso = length(enum.groups(gidx + 1).nodes(ni).srcs);
  for j = 1:nso
    enum.groups(gidx + 1).nodes(ni).srcs(j).optode_idx = enum.groups(gidx + 1).nodes(ni).srcs(j).optode_idx + 1;
  end
  
  ndo = length(enum.groups(gidx + 1).nodes(ni).dets);
  for j = 1:ndo
    enum.groups(gidx + 1).nodes(ni).dets(j).optode_idx = enum.groups(gidx + 1).nodes(ni).dets(j).optode_idx+ 1;
  end
 
  
end

% Remove extraneous non-canonical fields from the nodes structure
enum.groups(gidx + 1).nodes = rmfield(enum.groups(gidx + 1).nodes, ...
  {'busid', 'fwmajor', 'fwminor', 'fwpatch', 'idx', 'scqid', 'temp', 'type', 'vin'});

% Insert the source power into the nodal enumeration
%
% The canonical form does not accept souce powers which vary over time, so we must assert
% this to be true before we continue
for i = 1:size(srcpwr,3)
  
  if i == 1
    % srcpwr_ref is now [node x wavelength]
    srcpwr_ref = srcpwr(:,:,i);
  else
    if ~all(all(srcpwr_ref == srcpwr(:,:,i)))
      warning('Time varying source powers are not supported, output will be inconsistent');
    end
  end
  
end

% Insert as described
nd = length(enum.groups(gidx + 1).nodes);
for ni = 1:nd
  node = enum.groups(gidx +1).nodes(ni);
  nq = length(node.srcs);
  for qi = 1:nq
    src = enum.groups(gidx +1).nodes(ni).srcs(qi);
    
    % TODO: Wavelengths are hard-coded here, which is incorrect and should be asserted
    %       against at very least.
    if ni == 1
      if nq == 1
        warning('Assertion required on wavelengths for source power assignment');
      end
    end
    
    if(src.wl == 735)
      enum.groups(gidx +1).nodes(ni).srcs(qi).power = srcpwr_ref(ni, 1); 
    else
      enum.groups(gidx +1).nodes(ni).srcs(qi).power = srcpwr_ref(ni, 2);
    end
    
  end
  
end

%%% Check for bad frames
if(any(errcnt(:)))
    warning([...
        'The specified LUFR file contains incomplete frames. Consult the frame error '...
        'counter to determine those frames for which the data is incomplete. Note that '...
        'missing data will present as NaN in the associated data structure']);
end


%%% Add the layout
%
% LUFR files may contain embedded layouts, or we can take it from the user, or complain.
%


% Load layout or take from the user
if isempty(layout_override) && (n_layouts == 0)
  
  % The user has nor provided a layout, we must use the embedded data if it exists
  % enum.groups(gidx + 1).layout = [];
   warning([...
      'The specified LUFR file does not contain an embedded layout file, and no layout '...
      'has been specified when calling this function. The returned enumeration will '...
      'lack layout information, and it will not be possible to convert this file to '...
      'formats which require a layout. Specify an appropriate layout file to supress '... 
      'this warning.']);
      
else
    
  % The user has supplied a layout file
  if ischar(layout_override)
    
    try
      layout_raw = lumofile.read_layout(layout_override);
      enum.groups(gidx + 1).layout = layout_raw;
      fprintf('LUFR file using user-specified layout file\n');
    catch e
       fprintf('An error occurred loading the specified layout file %s\n', layout_override);
       rethrow e
    end
      
  elseif isstruct(layout_override)
    enum.groups(gidx + 1).layout = layout_override;
    fprintf('LUFR file using user-specified layout structure\n');
  elseif (n_layouts == 1)
    
    try
      layout_embed = jsondecode(layoutstr.');
    catch e
      fprintf('Error parsing embedded layout file\n');
      rethrow(e);
    end
    enum.groups(gidx + 1).layout = lumofile.proc_layout(layout_embed);
    
    fprintf('LUFR file using embedded layoutfile\n');
  else  
    error('The specified layout is neither a layout filename nor structure, consult help');
  end
    
end

if isfield(enum.groups(gidx + 1), 'layout')
  lumomat.validate_layout(enum, gidx + 1);
else
  enum.groups(gidx + 1).layout = [];
end

% Build event output structure 
ne = length(evtim);
if ne > 0
  for i = 1:ne
    % Note that events are stored in seconds
    events(i) = struct('mark', convertStringsToChars(evstr{i}), 'timestamp', evtim(i)*1e-3);
  end
else
  events = [];
end

% Filter events (remove non printable characters)
evfilt = true(ne,1);
for i = 1:ne
    if(length(events(i).mark) == 1)
        if events(i).mark < char(32)
            evfilt(i) = false;
        end
    end
end

events = events(evfilt);

if(sum(~evfilt) > 0)
    fprintf('Filtering %d/%d invalid event markers from recording\n', sum(~evfilt), length(evfilt));
end


% Build data output structure
data = struct('chn_dat', chdat, ...
              'chn_dt', tchdat*1e3, ...
              'chn_fps', fps, ...
              'chn_sat', satflag,...
              'err_cnt', errcnt,...
              'nframes', size(chdat, 2), ...
              'nchns', size(chdat, 1), ...
              'node_temp', tmpdat, ...
              'node_mpu_dt', 10e-3, ...
              'node_mpu_fps', 100, ...
              'node_acc', accdat, ...
              'node_gyr', gyrdat);
              
% Add the channel permutation
if(apply_filter)
  enum.groups(gidx+1).chn_perm = chperm;
  enum.groups(gidx+1).channels = enum.groups(gidx+1).channels(chperm);
end

% Remove all enumeration data other than that requested
enum.groups = enum.groups(gidx+1);

end
