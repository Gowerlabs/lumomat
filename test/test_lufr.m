% Test LUFR input
%
% Testing requires the lumofile test sample package, which is not public as it contains real
% data recordings from subjects.
%
%
%   (C) Gowerlabs Ltd., 2022
%

[cpath, ~, ~] = fileparts(mfilename('fullpath'));

lufr_sample_files = {...
  'sample_lufr_18_v2_1.lufr', ...
  'sample_lufr_17_v2_1.lufr'};

layout_sample_files = {...
  'layout_54_20155531.json',...
  'layout_54_20155531.json'};

%% Test 1: Convert from LumoData object to SNIRF

for i = 1:length(lufr_sample_files)

  % Get LUFR file and layout file
  lufr_sample_fn = fullfile(cpath, 'samples', lufr_sample_files{i});
  layout_sample_fn = fullfile(cpath, 'samples', layout_sample_files{i});

  % Load raw data
  ld = LumoData(lufr_sample_fn, 'layout', layout_sample_fn);

  [p,n,e] = fileparts(lufr_sample_fn);
  snirffn = fullfile(p, [n '.snirf']);

  snirf = ld.write_SNIRF(snirffn);

end

%% Test 2: Test big GID

[enum, data, events] = lumofile.read_lufr(fullfile(cpath, 'samples', 'sample_lufr_big_gid.lufr'));

%% Test 3: No events

[enum, data, events] = lumofile.read_lufr(fullfile(cpath, 'samples', 'sample_lufr_no_evt.lufr'));
ld = LumoData(fullfile(cpath, 'samples', 'sample_lufr_no_evt.lufr'), 'layout', fullfile(cpath, 'samples', 'layout_54_20155531.json'));

%% Test 4: Check sample against NIRS output
[enum, data, events] = lumofile.read_lufr(fullfile(cpath, 'samples', 'sample_lufr_short2.lufr'), ...
  'layout',  fullfile(cpath, 'samples', 'layout_54_20155531.json'));

nirs_test = lumofile.write_NIRS([], enum, data, events, 'sd_style', 'flat');

nirs_sample = load(fullfile(cpath, 'samples', 'sample_lufr_short2.nirs.mat'));

% Sample NIRS files contain an empty string and zero array when there are no events,
% whereas we have an empty cell and and empty variable. Test only when this is not the
% case
if(length(nirs_sample.CondNames) == 1)
  if(isempty(nirs_sample.CondNames{1}))
    test_s = false;
  end
else
  test_s = true;
end

% Basic variables
if test_s
  assert(all(all(nirs_sample.s == nirs_test.s)))
  assert(all([nirs_sample.CondNames{:}] == [nirs_test.CondNames{:}]))
end
assert(norm(nirs_sample.t - nirs_test.t) < 1e-9)

if norm(nirs_sample.d - nirs_test.d) > 1e-9
  warning('Data norm exceeds 1x10^-6');
end

% SD
assert(nirs_sample.SD.nSrcs == nirs_test.SD.nSrcs);
assert(nirs_sample.SD.nDets == nirs_test.SD.nDets);
assert(all(nirs_sample.SD.Lambda == nirs_test.SD.Lambda));
% assert(norm(nirs_sample.SD.SrcPos - nirs_test.SD.SrcPos) < 1e-9);
% assert(norm(nirs_sample.SD.DetPos - nirs_test.SD.DetPos) < 1e-9);
% assert(all(nirs_sample.SD.SrcPowers == nirs_test.SD.SrcPowers));
assert(all(all(nirs_sample.SD.MeasList == nirs_test.SD.MeasList)));
assert(all(nirs_sample.SD.MeasListAct == nirs_test.SD.MeasListAct));

% SD3D
assert(nirs_sample.SD3D.nSrcs == nirs_test.SD3D.nSrcs);
assert(nirs_sample.SD3D.nDets == nirs_test.SD3D.nDets);
assert(all(nirs_sample.SD3D.Lambda == nirs_test.SD3D.Lambda));
assert(norm(nirs_sample.SD3D.SrcPos - nirs_test.SD3D.SrcPos) < 1e-9);
assert(norm(nirs_sample.SD3D.DetPos - nirs_test.SD3D.DetPos) < 1e-9);
% assert(all(nirs_sample.SD3D.SrcPowers == nirs_test.SD3D.SrcPowers));
assert(all(all(nirs_sample.SD3D.MeasList == nirs_test.SD3D.MeasList)));
assert(all(nirs_sample.SD3D.MeasListAct == nirs_test.SD3D.MeasListAct));
% assert(all(nirs_sample.SD3D.MeasListActSat == nirs_test.SD3D.MeasListActSat));



%% Test 4: Load a version 3 as a LumoData object
lufr_ld = LumoData(fullfile(cpath, 'samples', 'sample_lufr_v3_with_layout.lufr'));