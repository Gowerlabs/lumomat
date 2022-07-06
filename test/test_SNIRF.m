% LumoData and lumofile SNIRF output tests
%
% Testing requires the lumofile test sample package, which is not public as it contains real
% data recordings from subjects.
%
%
%   (C) Gowerlabs Ltd., 2022
%

[path, ~, ~] = fileparts(mfilename('fullpath'));

lumo_sample_files = {...
  'sample_v011_1.LUMO',...
  'sample_v020_1.LUMO',...
  'sample_v030_1.LUMO',...
  'sample_v040_1.LUMO'};

nirs_sample_files = {...
  'sample_v011_1_flat.nirs.mat',...
  'sample_v020_1_flat.nirs.mat',...
  'sample_v030_1_flat.nirs.mat',...
  'sample_v040_1_flat.nirs.mat'};
  

lumo_sample_fn = fullfile(path, 'samples', lumo_sample_files);
nirs_sample_fn = fullfile(path, 'samples', nirs_sample_files);

layout_12_1_fn = fullfile(path, 'samples', 'layout_12_3735928559.json');

%% Test 1: Convert from LumoData object to SNIRF

for i = 1:length(lumo_sample_files)
  
  % Load raw data
  if i == 1
    % Special case for the first file, we insert a layout file
    ld = LumoData(lumo_sample_fn{1}, 'layout', layout_12_1_fn);
  else
    ld = LumoData(lumo_sample_fn{i});
  end
  
  [p,n,e] = fileparts(lumo_sample_fn{i});
  snirffn = fullfile(p, [n '.snirf']);
  
  ld.write_SNIRF(snirffn);
  
end

%%% Validation with PySNIRF2
%
% Recent versions of MATLAB can directly call Python in order to use the PySNIRF2 validation
% methods. To enable testing a suitable environment must be created, and activated within
% MATLAB. Additionally, the requirements (from requirements.txt) must be installed.
%
% For example, from the test directory:
%
% cd pyutil
% python3 -m venv ./pytenv
% ./pytenv/Scripts/Activate.ps1   # (Windows)
% pip install -r requirements.txt
%
% Then in MATLAB:
%
%   pyenv('Version', fullfile('.','pyutil','pytenv','Scripts','python.exe'))
% 
% Owing to limited MATLAB support this is not currently fleshed out, and instead you will
% have to manually validate the files, e.g., 
%
% >> from pysnirf2 import validateSnirf
% >> result = validateSnirf(r'samples\sample.snirf')
% >> assert result
