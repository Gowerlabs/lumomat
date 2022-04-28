function merge_layout(input_lumo_directory, layout_file, output_lumo_directory)
% LUMOFILE.MERGE_LAYOUT takes in a lumo file directory and a layout file directory, using
% them both to create a lumo file directory with an embedded layout file. If no output
% directory is defined then the given lumo_directory will be overwritten.
%
%   Parameters:
%   
%   input_lumo_directory    :   The path of the lumo directory to which the layout_file
%                               will be embedded into. WARNING: files in this folder will
%                               be overwirtten if no output_lumo_file_directory is
%                               defined.
%
%   layout_file             :   The path of the layout file which would be embedded into
%                               the lumo_directory, see `find_appdata_layout` to 
%                               obtain the layout file path.
%
%   output_lumo_directory   :   An optional parameter which sets the path to which the
%                               embedded lumo directory will be written to. If no output
%                               file is defined, then the input directory will be
%                               overwitten.
%

%% Check if output_lumo_directory exists, and if not, overwite

if ~exist('output_lumo_directory','var')
    warning('no output_layout_file_set, overwitting \"' + input_lumo_directory + '\"')
    output_lumo_directory = input_lumo_directory;
end

%% Check if the 1st 2 directories exist and the last one doesn't.

if(~exist(input_lumo_directory,'dir'))
    error('LUMO file (%s): invalid directory not found\n', input_lumo_directory);
end

if(exist(output_lumo_directory,'dir') && output_lumo_directory ~= input_lumo_directory)
    error('LUMO file (%s): directory already exists\n', input_lumo_directory);
end

if(~isfile(layout_file))
    error('Layout file (%s): invalid directory not found\n', layout_file);
end

%% Reading and opening the Meta Data file

metadata_fn = fullfile(input_lumo_directory, 'metadata.toml');

if exist(metadata_fn, 'file') ~= 2
  error('LUMO file (%s): invalid metadata not found\n', input_lumo_directory);
end

try
  raw = fileread(metadata_fn);
  metadata = lumofile.toml.decode(raw);
catch e
  fprintf(2, 'LUMO file (%s): error parsing metadata file %s\n', input_lumo_directory, metadata_fn);
  rethrow(e);
end

output_metadata = metadata;

%% Parse version information

lf_missing_layout_ver = [0 0 1; 0 1 0; 0 1 1; 0 4 0];

lf_ver = reqfield(metadata, 'lumo_file_version', input_lumo_directory);
try
  lf_ver_num = str2double(strsplit(lf_ver, '.'));
catch e
  fprintf(2, 'LUMO file (%s): error parsing fiile version number\n', input_lumo_directory);
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
    warning('This lumo file version should already contain a layout file.');
end

%% Set output file layout_file to hardware_file if applicable.

lf_meta_fns = reqfield(metadata, 'file_names');

if ~isfield(lf_meta_fns, 'hardware_file')
    output_metadata.file_names.hardware_file = lf_meta_fns.lf_meta_fns;
end
