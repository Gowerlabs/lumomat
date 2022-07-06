function lumo_merged_fn = merge_layout(lumo_in_fn, layout_file)
% LUMOFILE.MERGE_LAYOUT Merge a layout file with an existing LUMO recording
%
% lumo_merged_filename = LUMOFILE.MERGE_LAYOUT(lumo_filename, layout_filename)
%
% LUMOFILE.MERGE_LAYOUT copies the input LUMO file and merges the specified layout file into
% the copy. The merged output file name is returned to the user.
%
%
%   Parameters:
%   
%   lumo_filename:  Filename of the LUMO file into which the layout file will be embedded.
%
%   layout_file:    Filename of the layout file to be merged. If no filename is provided the
%                   function will search for the layout file and print its path.
%
% Notes:
%
% The purpose of this function is to merge an appropriate layout file into a LUMO recording
% which lacks a layout. This may occurr since some versions of the LUMOview software do not
% automically embed the layout file. 
%
% If the intent is instead to merge a subject specific layout file, this is better achieved
% by providing the alternative layout file when loading the data.
%
%
%   (C) Gowerlabs Ltd., 2022
%

% Normalise strings
[lumo_in_fn, layout_file] = convertStringsToChars(lumo_in_fn, layout_file);

%%% Check for existance of the input LUMO directory and layout file
if exist(lumo_in_fn, 'dir') ~= 7
  error('The specified LUMO file (%s) cannot be found', lumo_in_fn)
else
  fprintf('Found LUMO file %s\n', lumo_in_fn);
end

if nargin > 1

  if exist(layout_file, 'file') ~= 2
    error('The specified layout file (%s) cannot be found', layout_file)
  else
    fprintf('Found layout file %s\n', layout_file);
  end
  layout_specced = true;
else
  fprintf('No layout filename specified, will search for layout file in default locations\n');
  layout_specced = false;
end
%%% Get the input metadata
metadata_fn = fullfile(lumo_in_fn, 'metadata.toml');

if exist(metadata_fn, 'file') ~= 2
  error('LUMO file (%s): invalid metadata not found\n', lumo_in_fn);
end

try
  raw = fileread(metadata_fn);
  metadata = lumofile.toml.decode(raw);
catch e
  fprintf(2, 'LUMO file (%s): error parsing metadata file %s\n', lumo_in_fn, metadata_fn);
  rethrow(e);
end

output_metadata = metadata;

%%% Parse and check version against supported
lf_known_ver = [0 0 1; 0 1 0; 0 1 1; 0 2 0; 0 3 0; 0 4 0; 0 5 0];
lf_mergable_ver = [0 2 0; 0 3 0; 0 4 0; 0 5 0];

lf_ver = reqfield(metadata, 'lumo_file_version', lumo_in_fn);
try
  lf_ver_num = str2double(strsplit(lf_ver, '.'));
catch e
  fprintf(2, 'LUMO file (%s): error parsing fiile version number\n', lumo_in_fn);
  rethrow(e);
end

% We need to get the (required) file name field here, in order to deal with a version number
% ambiguity whereby 0.2.0 files can report as 0.1.1.
lf_meta_fns = reqfield(metadata, 'file_names');
if all(lf_ver_num == [0 1 1])
  if isfield(lf_meta_fns, 'hardware_file')
    lf_ver_num = [0 2 0];
  end
end

if ~ismember(lf_known_ver, lf_ver_num, 'rows')
  error('LUMO file (%s): version %d.%d.%d is not supported by this software version\n', ...
    lumo_in_fn, lf_ver_num(1), lf_ver_num(2), lf_ver_num(3));
else
  fprintf('LUMO file version %d.%d.%d\n', ...
    lf_ver_num(1), lf_ver_num(2), lf_ver_num(3));
end

% If we were to process version v0.0.1 0 v0.1.1 we would need to reformat the metadata to
% allow the layout filename to be stored in the appropriate field.
%
% The only way to do this is to upgrade the file to v0.2.0 or above, which this script
% currently does not do. Instead, we error out in this circumstance.
%
% % Set output file layout_file to hardware_file if applicable.
% lf_meta_fns = reqfield(metadata, 'file_names');
% 
% if ~isfield(lf_meta_fns, 'hardware_file')
%     output_metadata.file_names.hardware_file = lf_meta_fns.layout_file;
% end

if ~ismember(lf_mergable_ver, lf_ver_num, 'rows')
  error('Unable to merge layouts into LUMO files of version, contact Gowerlabs');
end


%%% Either search for the layour file or load and continue
[lf_desc] = load_lumo_desc(lumo_in_fn);
[enum, ~] = load_lumo_enum(lumo_in_fn, lf_desc);

if length(enum.groups) > 1
  error(['Input LUMO file (%s) contains recordings from more than one group, unable to'...
         'continue']);
end

[uid_hex, uid_name] = lumomat.norm_gid(enum.groups(1).id);

if ~layout_specced
  layout_installed_fn = lumofile.find_install_layout(hex2dec(uid_hex));
  if isempty(layout_installed_fn)
    error(['No layout file was specified for merging, and a layout file cannot be ' ...
           'automatically found for this LUMO file (group %s / %s)\n'], ...
           uid_hex, uid_name);
  else
    error(['\nNo layout file was specified for merging, but a suitable layout file has '...
             'been located in the standard installation directory: %s\n'], ...
             layout_installed_fn);
  end
end

layout_data = lumofile.read_layout(layout_file);

if ~strcmp(uid_hex, layout_data.id)
  warning(['Specified layout (group %s / %s) does not match LUMO file (group %s /% s), '...
           'output file may be inconsistent'], ...
           layout_data.id, layout_data.name, uid_hex, uid_name);
end


%%% Create output directory and copy contents
[lf_filepath, lf_name, lf_ext] = fileparts(lumo_in_fn);
lumo_merged_fn = fullfile(lf_filepath, [lf_name '_merged_' uid_name lf_ext]);
fprintf('Writing merged file to %s\n', lumo_merged_fn);

[status,msg] = copyfile(lumo_in_fn, lumo_merged_fn);

if ~status
  fprintf(2, 'Error duplicating the input LUMO file: %s\n', msg);
  return;
end

%%% (over)write layout.json
[status, msg] = copyfile(layout_file, fullfile(lumo_merged_fn, 'layout.json'));
if(~status)
    error('Error merging layout file (%s) into output LUMO file (%s): %s', ...
             layout_file, lumo_merged_fn, msg);
end

% Set and overwrite metadat
output_metadata.file_names.layout_file = 'layout.json';

try  
    fid = fopen(fullfile(lumo_merged_fn, "metadata.toml"), 'w');
    fprintf(fid, '%s', lumofile.toml.encode(output_metadata));
    fclose(fid);
catch e
    fprintf('Error updating metadata in output LUMO file (%s)\n', lumo_merged_fn);
    rethrow(e);
end
    

fprintf('Layout merge complete\n');

end

