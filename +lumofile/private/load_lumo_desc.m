function [lf_desc] = load_lumo_desc(lf_dir)
% LOAD_LUMO_DESC
%
% Create the lumo file description, returning a structure of information required for
% parsing, including the version, file-names, and flags indicating the available contents.
% All files referenced in the structure have been confirmed to exist within the file.
%

% 0. Some constants
lf_known_ver = [0 0 1; 0 1 0; 0 1 1; 0 2 0; 0 3 0; 0 4 0; 0 5 0];

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
