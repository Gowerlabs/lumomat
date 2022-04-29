% Layout location find and merging tests
%
% Check that layout files can be located when available, and that merged output remains
% identical to original files.%
%
%
%   (C) Gowerlabs Ltd., 2022
%

[lmpath, ~, ~] = fileparts(mfilename('fullpath'));

% Cross-check datasets
lumo_xc_sample_files = {...
  'sample_v011_1.LUMO',...
  'sample_v020_1.LUMO',...
  'sample_v030_1.LUMO',...
  'sample_v040_1.LUMO'};

nirs_xc_sample_files = {...
  'sample_v011_1_flat.nirs.mat',...
  'sample_v020_1_flat.nirs.mat',...
  'sample_v030_1_flat.nirs.mat',...
  'sample_v040_1_flat.nirs.mat'};

% Independent version samples
lumo_idpver_samples = {...
  'sample_v050_1.LUMO'};

lumo_xc_sample_fn = fullfile(lmpath, 'samples', lumo_xc_sample_files);
nirs_xc_sample_fn = fullfile(lmpath, 'samples', nirs_xc_sample_files);

lumo_idpver_sample_fn = fullfile(lmpath, 'samples', lumo_idpver_samples);


layout_12_1_fn = fullfile(lmpath, 'samples', 'layout_12_3735928559.json');


%% Test 1: Map forwards and backwards

group_id_samples = [742390278; 20155499; 12];

for i = 1:length(group_id_samples)
  
  [hex, name] = lumomat.norm_gid(group_id_samples(i));
  assert(group_id_samples(i) == lumomat.norm_gname(name))
 
end

%% Test 2: check for layout present warning

for i = 1:length(lumo_xc_sample_fn)
  [enum, data,  events] = lumofile.read_lumo(lumo_xc_sample_fn{1});
end


%% Test 3: merge layout files
try
  lumofile.merge_layout(fullfile(lmpath, 'samples', 'sample_v040_no_layout.LUMO'));
  error('Merge layout should throw');
catch e
end
  



lumo_merged_fn = lumofile.merge_layout(fullfile(lmpath, 'samples', 'sample_v040_no_layout.LUMO'), ...
                                       'C:\Users\Sam\AppData\Local\Gowerlabs\Lumo\coordinates_268.json');
ld = LumoData(lumo_merged_fn);

assert(isstruct(ld.enum.groups.layout))




