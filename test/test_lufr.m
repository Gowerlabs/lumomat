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
  'sample_lufr_18_v2_1.lufr', ...
  'sample_lufr_17_v2_1.lufr'};

layout_sample_files = {...
  'layout_54_20155531.json',...
  'layout_54_20155531.json'};

%% Test 1: Convert from LumoData object to SNIRF

for i = 1:length(lufr_sample_files)

  % Get LUFR file and layout file
  lufr_sample_fn = fullfile(path, 'samples', lufr_sample_files{i});
  layout_sample_fn = fullfile(path, 'samples', layout_sample_files{i});
  
  % Load raw data
  ld = LumoData(lufr_sample_fn, 'layout', layout_sample_fn);

  [p,n,e] = fileparts(lufr_sample_fn);
  snirffn = fullfile(p, [n '.snirf']);
  
  snirf = ld.write_SNIRF(snirffn, 'style', 'mne-nirs');
  
end

%% Test 2: Test big GID

[enum, data, events] = lumofile.read_lufr(fullfile(path, 'samples', 'sample_lufr_big_gid.lufr'));

%% Test 3: No events

[enum, data, events] = lumofile.read_lufr(fullfile(path, 'samples', 'sample_lufr_no_evt.lufr'));
ld = LumoData(fullfile(path, 'samples', 'sample_lufr_no_evt.lufr'), 'layout', fullfile(path, 'samples', 'layout_54_20155531.json'));