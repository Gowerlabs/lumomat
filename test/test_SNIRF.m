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

layout_12_1_fn = fullfile(path, 'samples', 'layout_12_1.json');

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

%%% TODO
%
% - Automatic SNIRF validation
%
% Manual validation (windows):
%  .\env\Scripts\activate.ps1
%
% >> from pysnirf2 import validateSnirf
% >> result = validateSnirf(r'samples\sample.snirf.snirf')
% >> assert result








