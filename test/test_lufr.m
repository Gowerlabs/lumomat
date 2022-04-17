% LumoData and lumofile LUFR input
%
% Testing requires the lumofile test sample package, which is not public as it contains real
% data recordings from subjects.
%
%
%   (C) Gowerlabs Ltd., 2022
%

[path, ~, ~] = fileparts(mfilename('fullpath'));

lufr_sample_files = {...
  'sample_lufr_v200_1.lufr'};

lumo_sample_fn = fullfile(path, 'samples', lufr_sample_files);

laout_sample_v200_1 = fullfile(path, 'samples', 'layout_sample_lufr_v200_1.json');

%% Test 1: Convert from LumoData object to SNIRF

for i = 1:length(lufr_sample_files)
  
  % Load raw data
  ld = LumoData(lumo_sample_fn{1}, 'layout', laout_sample_v200_1);

  [p,n,e] = fileparts(lumo_sample_fn{i});
  snirffn = fullfile(p, [n '.snirf']);
  
  snirf = ld.write_SNIRF(snirffn, 'ordering', 'mne-nirs');
  
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
% >> result = validateSnirf(r'samples\sample.snirf.snirf')
% >> assert result
