% lumofile tests
%
% Testing requires the lumofile test sample package, which is not public as it contains real
% data recordings from subjects.
%
% Read tests:
%
% - Lumo file not found
% - No metadata in the file
% - Missing intensity files (referenced)
% - No intensity files
% - Garbage metadata file
% - Garbage hardware file
% - Garbange events
% - Garbage raw data, incl.
%   - version mismatch
%   - endienness mismatch
%   - mismatch between reported size
% - Checks on enumeration mapping
% - Checks on data location mapping
% - Parse checks on v0.0.1, v0.1.0, v0.2.0, v0.3.0, v0.4.0
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


%% Test 1: parse all supported lumo file versions

for i = 1:length(lumo_sample_files)
  [enum, data,  events] = lumofile.read_lumo(lumo_sample_fn{i});
end

%% Test 2b: Check v0.2.0 upwards against NIRS with embedded JSON

for i = 1:length(lumo_sample_files)

  nirs_sample = load(nirs_sample_fn{i});

  if i == 1
    % Special case for the first file, we insert a layout file
    [enum, data, events] = lumofile.read_lumo(lumo_sample_fn{1}, 'layout', layout_12_1_fn);
  else
    [enum, data, events] = lumofile.read_lumo(lumo_sample_fn{i});
  end
  
  % Build the 'flat' layout and test
  nirs_test = lumofile.write_NIRS([], enum, data, events, 'sdstyle', 'flat');

  % Basic variables
  assert(all(nirs_sample.s == nirs_test.s))
  assert(all(nirs_sample.t == nirs_test.t))
  assert(all(nirs_sample.d == nirs_test.d))
  assert(all(nirs_sample.CondNames{:} == nirs_test.CondNames{:}))

  % SD
  assert(nirs_sample.SD.nSrcs == nirs_test.SD.nSrcs);
  assert(nirs_sample.SD.nDets == nirs_test.SD.nDets);
  assert(all(nirs_sample.SD.Lambda == nirs_test.SD.Lambda));
  assert(all(nirs_sample.SD.SrcPos == nirs_test.SD.SrcPos));
  assert(all(nirs_sample.SD.DetPos == nirs_test.SD.DetPos));
  assert(all(nirs_sample.SD.SrcPowers == nirs_test.SD.SrcPowers));
  assert(all(nirs_sample.SD.MeasList == nirs_test.SD.MeasList));
  assert(all(nirs_sample.SD.MeasListAct == nirs_test.SD.MeasListAct));
  
  % SD3D
  assert(nirs_sample.SD3D.nSrcs == nirs_test.SD3D.nSrcs);
  assert(nirs_sample.SD3D.nDets == nirs_test.SD3D.nDets);
  assert(all(nirs_sample.SD3D.Lambda == nirs_test.SD3D.Lambda));
  assert(all(nirs_sample.SD3D.SrcPos == nirs_test.SD3D.SrcPos));
  assert(all(nirs_sample.SD3D.DetPos == nirs_test.SD3D.DetPos));
  assert(all(nirs_sample.SD3D.SrcPowers == nirs_test.SD3D.SrcPowers));
  assert(all(nirs_sample.SD3D.MeasList == nirs_test.SD3D.MeasList));
  assert(all(nirs_sample.SD3D.MeasListAct == nirs_test.SD3D.MeasListAct));
  
  % Build the 'standard' layout and test
  nirs_test = lumofile.write_NIRS([], enum, data, events);

  % Basic variables
  assert(all(nirs_sample.s == nirs_test.s))
  assert(all(nirs_sample.t == nirs_test.t))
  assert(all(nirs_sample.d == nirs_test.d))
  assert(all(nirs_sample.CondNames{:} == nirs_test.CondNames{:}))

  % SD - note that we can compare the sample SD3D to the test SD here because the standard
  % format uses the 3D in th main structure.
  assert(nirs_sample.SD3D.nSrcs == nirs_test.SD.nSrcs);
  assert(nirs_sample.SD3D.nDets == nirs_test.SD.nDets);
  assert(all(nirs_sample.SD3D.Lambda == nirs_test.SD.Lambda));
  assert(all(nirs_sample.SD3D.SrcPos == nirs_test.SD.SrcPos));
  assert(all(nirs_sample.SD3D.DetPos == nirs_test.SD.DetPos));
  assert(all(nirs_sample.SD3D.SrcPowers == nirs_test.SD.SrcPowers));
  assert(all(nirs_sample.SD3D.MeasList == nirs_test.SD.MeasList));
  assert(all(nirs_sample.SD3D.MeasListAct == nirs_test.SD.MeasListAct));
  
end




