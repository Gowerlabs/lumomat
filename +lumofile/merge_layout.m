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

%% Reading and opening the Meta Data file

metadata_fn = fullfile(input_lumo_directory, 'metadata.toml');

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

output_metadata = metadata;

