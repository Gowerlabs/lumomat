% LumoData and lumofile tests
%
% Testing requires the lumofile test sample package, which is not public as it contains real
% data recordings from subjects.
%
%
%   (C) Gowerlabs Ltd., 2022
%

%%% TODO
%
% Additonal tests required
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


%% Test 1: parse all supported lumo file versions

for i = 1:length(lumo_xc_sample_fn)
  [enum, data,  events] = lumofile.read_lumo(lumo_xc_sample_fn{i});
end

%% Test 2: Check against NIRS

for i = 1:length(lumo_xc_sample_files)

  nirs_sample = load(nirs_xc_sample_fn{i});

  if i == 1
    
    % Special case for the first file, we insert a layout file
    [enum, data, events] = lumofile.read_lumo(lumo_xc_sample_fn{1}, 'layout', layout_12_1_fn);
  else
    [enum, data, events] = lumofile.read_lumo(lumo_xc_sample_fn{i});
  end
  
  % Build the 'flat' layout and test
  nirs_test = lumofile.write_NIRS([], enum, data, events, 'sd_style', 'flat');
  
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
    assert(all(nirs_sample.s == nirs_test.s))
    assert(all(nirs_sample.CondNames{:} == nirs_test.CondNames{:}))
  end
  assert(norm(nirs_sample.t - nirs_test.t) < 1e-9)
  
  if norm(nirs_sample.d - nirs_test.d) > 1e-6
    warning('Data norm exceeds 1x10^-6');
  end
  
  assert(norm(nirs_sample.d - nirs_test.d) < 1)
    
  % SD
  assert(nirs_sample.SD.nSrcs == nirs_test.SD.nSrcs);
  assert(nirs_sample.SD.nDets == nirs_test.SD.nDets);
  assert(all(nirs_sample.SD.Lambda == nirs_test.SD.Lambda));
  assert(norm(nirs_sample.SD.SrcPos - nirs_test.SD.SrcPos) < 1e-9);
  assert(norm(nirs_sample.SD.DetPos - nirs_test.SD.DetPos) < 1e-9);
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
  assert(all(nirs_sample.SD3D.MeasListActSat == nirs_test.SD3D.MeasListActSat));
  
  % Build the 'standard' layout and test
  nirs_test = lumofile.write_NIRS([], enum, data, events);

  % Basic variables
  if test_s
    assert(all(nirs_sample.s == nirs_test.s))
    assert(all(nirs_sample.CondNames{:} == nirs_test.CondNames{:}))
  end
  assert(norm(nirs_sample.t - nirs_test.t) < 1e-9)
  
  if norm(nirs_sample.d - nirs_test.d) > 1e-6
    warning('Data norm exceeds 1x10^-6');
  end
  
  % SD - note that we can compare the sample SD3D to the test SD here because the standard
  % format uses the 3D in th main structure.
  assert(nirs_sample.SD3D.nSrcs == nirs_test.SD.nSrcs);
  assert(nirs_sample.SD3D.nDets == nirs_test.SD.nDets);
  assert(all(nirs_sample.SD3D.Lambda == nirs_test.SD.Lambda));
  assert(all(all(nirs_sample.SD3D.SrcPos == nirs_test.SD.SrcPos)));
  assert(all(all(nirs_sample.SD3D.DetPos == nirs_test.SD.DetPos)));
  % assert(all(nirs_sample.SD3D.SrcPowers == nirs_test.SD.SrcPowers));
  assert(all(all(nirs_sample.SD3D.MeasList == nirs_test.SD.MeasList)));
  assert(all(nirs_sample.SD3D.MeasListAct == nirs_test.SD.MeasListAct));
  assert(all(nirs_sample.SD3D.MeasListActSat == nirs_test.SD.MeasListActSat));
  
end

%% Test 3: Test high level interface

for i = 1:length(lumo_xc_sample_files)
  
  % Load raw data
  if i == 1
    % Special case for the first file, we insert a layout file
    ld = LumoData(lumo_xc_sample_fn{1}, 'layout', layout_12_1_fn);
  else
    ld = LumoData(lumo_xc_sample_fn{i});
  end
  
  % Load the NIRS sample
  nirs_sample = load(nirs_xc_sample_fn{i});
  
  % Convert to NIRS to get sorting permutation
  nirs_test = lumofile.write_NIRS([], ld.enum, ld.data, ld.evts, 'sd_style', 'flat');
  
    % Check vectors such as time
  ml_t = nirs_sample.t; 
  ld_t = ld.chn_time;
  assert(norm(ml_t - ld_t.') < 1e-9);
  
      
  for cidx = 1:length(ld.enum.groups.channels)
    
    % Get canonical channel indices
    lidx = nirs_test.lumo.chn_sort_perm(cidx);
    
    % Get the channel information
    ld_info = ld.chn_info(lidx);
    
    % Get NIRS parameters via the MeasList
    ml = nirs_sample.SD3D.MeasList(cidx,:);        
    ml_srcpos_3D = nirs_sample.SD3D.SrcPos(ml(1),:);
    ml_detpos_3D = nirs_sample.SD3D.DetPos(ml(2),:);
    ml_srcpos_2D = nirs_sample.SD.SrcPos(ml(1),:);
    ml_detpos_2D = nirs_sample.SD.DetPos(ml(2),:);    
    
    % Assert
    assert(all(ld_info.src_coords_3d == ml_srcpos_3D));
    assert(all(ld_info.det_coords_3d == ml_detpos_3D));
    assert(all([ld_info.src_coords_2d 0] == ml_srcpos_2D));
    assert(all([ld_info.det_coords_2d 0] == ml_detpos_2D));
  
    % Manually find some channels
    % Search for some channels and check these match the enumeration
    % Construct a global index and check against the NIRS output
    
  end
  
end

%% Test 4: New versions load
for i = 1:length(lumo_idpver_sample_fn)
  [enum, data,  events] = lumofile.read_lumo(lumo_idpver_sample_fn{i});
end

%% Test 5: Bad layout merging

%% Test 6: Non spectroscopic

  





