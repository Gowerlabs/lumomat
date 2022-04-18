% Test LUFR input
%
% Testing requires the lumofile test sample package, which is not public as it contains real
% data recordings from subjects.
%
%
%   (C) Gowerlabs Ltd., 2022
%

[path, ~, ~] = fileparts(mfilename('fullpath'));

lufr_sample_files = {...
  'sample_lufr_v2_1.lufr', ...
  'sample_lufr_18_v2_1.lufr'};

lumo_sample_fn = fullfile(path, 'samples', lufr_sample_files);

layout_fn = {fullfile(path, 'samples', 'layout_sample_lufr_v2_1.json'), ...
  fullfile(path, 'samples', 'layout_20155531.json')};

%% Test 1: Convert from LumoData object to SNIRF

for i = 1:length(lufr_sample_files)
  
  % Load raw data
  ld = LumoData(lumo_sample_fn{i}, 'layout', layout_fn{i});

  [p,n,e] = fileparts(lumo_sample_fn{i});
  snirffn = fullfile(p, [n '.snirf']);
  
  snirf = ld.write_SNIRF(snirffn, 'style', 'mne-nirs');
  
end

