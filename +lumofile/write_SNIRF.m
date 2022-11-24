function [snirf] = write_SNIRF(snirffn, enum, data, events, varargin)
% LUMOFILE.WRITE_SNIRF Write LUMO data to SNIRF format
%
%   [snirf] = LUMOFILE.WRITE_SNIRF(snirffn, enum, data, events)
%
%   LUMOFILE.WRITE_SNIRF writes to LUMO data to disk in the SNIRF v1.0 format, according
%   to the specification:
%
%   https://github.com/fNIRS/snirf/blob/v1.0/snirf_specification.md
%
%   Paramters:
%
%     snirffn:              The file name of the output SNIRF file.
%
%     enum, data, events:   Data structures returned by LUMOFILE.READ
%
%   Optional Parameters:
%
%   'meta':     Additional metadata can be written to the SNIRF file containing an 
%               abbreviated and flattened representation of the canonical layout for each
%               group defined in `/nirs(i)`.
%
%               Enabling this option will result in a warning from the official SNIRF 
%               validation utility.
%
%               'standard':   (default) no extended metadata   
%               'extended':   write additional metadata for `/nirs(i)` in `/lumo(i)`
%
%   Returns:
%
%     snirf:    A representaiton of the SNIRF file in structure format.
%
%   Details:
%
%   The SNIRF output is described by specification. Certain addition lumo specific metadata
%   is included, details of which can be found in the documentation for this package.
%
%   To construct a SNIRF output description of the system, the LUMO enumeration is
%   transformed from the canonical node local format to a globally indexed spectroscopic
%   format, such that channels are described by the tuple:
%
%   (src_position(i), det_postion(j), src_wavelength(k))
%
%   The canonical LUMO format permits more complex representations of the system than can be
%   accommodated by this format, so the actual enumeration is checked to ensure that it can
%   conform to this representation.
%
%   When constructing the global enumeration, the template layout (e.g. collection of all
%   docks in a group/cap) is collapsed to form a description containing only the occupied
%   docks. Extended information is provided containing the complete dock layout.
%
%
% See also LUMO_READ
%
%
%   (C) Gowerlabs Ltd., 2022
%

%%% TODO
%
% Pending: https://github.com/fNIRS/snirf/issues/103, implement modified channel descriptors
%          for significant file size reduction.
%

% Normalise strings
snirffn = convertStringsToChars(snirffn);

% Get optional inputs
p = inputParser;
expected_styles_meta = {'standard', 'extended'};
addOptional(p, 'meta', 'standard', @(x) any(validatestring(x, expected_styles_meta)));
parse(p, varargin{:})
meta_style = p.Results.meta;


ng = length(enum.groups);

% Check that a layout is available for every group
for gidx = 1:ng
  if isempty(enum.groups(gidx).layout)
    error('Input LUMO enumeration does not contain an embedded layout file');
  end
end

fprintf('Writing SNIRF file %s...\n', snirffn);

% Open the file
try
  fid = H5F.create(snirffn, 'H5F_ACC_TRUNC', H5P.create('H5P_FILE_CREATE'), H5P.create('H5P_FILE_ACCESS'));
catch e
  fprintf('Error creating SNIRF file %s\n', snirffn);
  rethrow(e);
end

% Write format version
%
snirf.formatVersion = write_var_string(fid, '/formatVersion', '1.1');

% Over each group
%
for gidx = 1:ng
  
  % Build the global spectroscopic mapping
  %
  try
    [glch, glsrc, gldet, glwl] = lumofile.map_gs(enum, 'group', gidx);
  catch e
    fprintf('Error forming gloabl spectroscopic mapping\n');
    rethrow(e)
  end
  
  % Build permutation for ordering
  chn_perm = [];
    
  % Create NIRS root
  % /nirs{i}
  %
  if ng > 1
    nirs_group = create_group(fid, ['nirs' num2str(gidx)]);
  else
    nirs_group = create_group(fid, ['nirs']);
  end
  
  % Create metadata
  % /nirs{i}/metaDataTags
  %
  
  % Required fields
  nirs_meta_group = create_group(nirs_group, 'metaDataTags');
  mdt.SubjectID = write_var_string(nirs_meta_group, 'SubjectID', 'Subject Unknown');
  mdt.MeasurementDate = write_var_string(nirs_meta_group, 'MeasurementDate', 'unknown');
  mdt.MeasurementTime = write_var_string(nirs_meta_group, 'MeasurementTime', 'unknown');
  mdt.LengthUnit = write_var_string(nirs_meta_group, 'LengthUnit', 'mm');
  mdt.TimeUnit = write_var_string(nirs_meta_group, 'TimeUnit', 'ms');
  mdt.FrequencyUnit = write_var_string(nirs_meta_group, 'FrequencyUnit', 'Hz');
  
  % Optional fields
  mdt.sourcePowerUnit = write_var_string(nirs_meta_group, 'sourcePowerUnit', 'percent');
  
  % LUMO sepcific fields
  mdt.ManufacturerName = write_var_string(nirs_meta_group, 'ManufacturerName', 'Gowerlabs');
  mdt.Model = write_var_string(nirs_meta_group, 'Model', 'LUMO');  
  mdt.lumomatVersion = write_var_string(nirs_meta_group, 'lumomatVersion', lumomat.ver());
  mdt.groupID = write_var_string(nirs_meta_group, 'groupID', enum.groups(gidx).id);
  mdt.groupName = write_var_string(nirs_meta_group, 'groupName', enum.groups(gidx).name); 
  
  % Write global saturation
  if isfield(data(gidx), 'chn_sat')
    mdt.saturationFlags = write_int32(nirs_meta_group, 'saturationFlags', any(data(gidx).chn_sat, 2));
  end
  
  % Write frame errors
  if isfield(data(gidx), 'err_cnt')
    mdt.errorFlags = write_int32(nirs_meta_group, 'errorFlags',  logical(data(gidx).err_cnt > 0));
  end
    
  if strcmp(meta_style, 'extended')
    
    % LUMO specific metadata
    if ng > 1
      lumo_md_group = create_group(fid, ['lumo' num2str(gidx)]);
    else
      lumo_md_group = create_group(fid, ['lumo']);
    end
  
    %%% Output the canonical map
    dockmap = enum.groups(gidx).layout.dockmap;
    nc = length(enum.groups(gidx).channels);
    canmap = zeros(nc, 7, 'int32');
    canmap(:,1) = [enum.groups(gidx).channels.src_node_idx];  % Source node index
    canmap(:,2) = dockmap(canmap(:,1));                       % Source dock index

    src_idx = [enum.groups(gidx).channels.src_idx];
    for ci = 1:nc 
      src = enum.groups(gidx).nodes(canmap(ci,1)).srcs(src_idx(ci));
      canmap(:,3) = src.optode_idx;   % Optode index
      canmap(:,4) = src.wl;           % Wavelength
    end

    canmap(:,5) = [enum.groups(gidx).channels.det_node_idx];  % Detecot node index
    canmap(:,6) = dockmap(canmap(:,5));                       % Detector dock index

    det_idx = [enum.groups(gidx).channels.det_idx];
    for ci = 1:nc 
      det = enum.groups(gidx).nodes(canmap(ci,1)).dets(det_idx(ci));
      canmap(:,7) = det.optode_idx;   % Optode index
    end

    if ~isempty(chn_perm)
      canmap = canmap(chn_perm,:);
    end

    mdt.lumo.canonicalMap = write_int32(lumo_md_group, 'canonincalMap', canmap);


    %% Output abbreviated nodal enumeration
    nn = length(enum.groups(gidx).nodes);
    for i = 1:nn
      node_group = create_group(lumo_md_group, ['node' num2str(i)]);
      mdt.lumo.nodes(i).id = write_int32(node_group, 'id', enum.groups(gidx).nodes(i).id);
      mdt.lumo.nodes(i).revision = write_int32(node_group, 'revision', enum.groups(gidx).nodes(i).rev);
      mdt.lumo.nodes(i).firmwareVersion = write_var_string(node_group, 'firmwareVersion',  enum.groups(gidx).nodes(i).fwver);    
      H5G.close(node_group); 
    end

    %% Output abbreviated dock information
    nd = length(enum.groups(gidx).layout.docks);
    for i = 1:nd
      dock = enum.groups(gidx).layout.docks(i);
      dock_group = create_group(lumo_md_group, ['dock' num2str(i)]);
      mdt.lumo.docks(i).id = write_int32(dock_group, 'id', dock.id);
      mdt.lumo.docks(i).optodeNames = write_var_string(dock_group, 'optodeNames', {dock.optodes.name});

      % Write out the optode positions
      no = length(dock.optodes);
      optodePos2D = zeros(no, 2);
      optodePos3D = zeros(no, 3);
      for j = 1:no
        optode = dock.optodes(j);
        optodePos2D(j,:) = [optode.coords_2d.x optode.coords_2d.y];
        optodePos3D(j,:) = [optode.coords_3d.x optode.coords_3d.y optode.coords_3d.z];
      end

      mdt.lumo.docks(i).optodePos2D = write_double(dock_group, 'optodePos2D', optodePos2D);
      mdt.lumo.docks(i).optodePos3D = write_double(dock_group, 'optodePos3D', optodePos3D);  

      H5G.close(dock_group);
    end

  

    H5G.close(lumo_md_group);  
  end
  
  % Assign metadata tags
  snirf.nirs(gidx).metaDataTags = mdt;
    
  H5G.close(nirs_meta_group);
  
  %% Create data block
  % /nirs{i}/data{i}
  %
  nirs_data_group = create_group(nirs_group, 'data1');                                            % Note: single data block
  dt.dataTimeSeries = write_chn_dat_block(nirs_data_group, chn_perm, data(gidx).chn_dat); % Write data /nirs{i}/data1
  dt.time = write_double(nirs_data_group, 'time', [0 data(gidx).chn_dt]);                 % Write /nirs{i}/time
  dt.measurementList = write_measlist(nirs_data_group, chn_perm, gidx, enum, glch);       % Write measurementList 
  
  % Assign data block
  snirf.nirs(gidx).data(1) = dt;
  H5G.close(nirs_data_group);
  
  %% Create probe block
  % /nirs{i}/probe
  %
  nirs_probe_group = create_group(nirs_group, 'probe');
  snirf.nirs(gidx).probe = write_probe(nirs_probe_group, enum, gidx, glsrc, gldet, glwl);
  H5G.close(nirs_probe_group);
  
  %% Create stimulus group
  %
  if ~isempty(events)
    
    % Find all of the conditions and their indices    
    [markers, ~, marker_perm] = unique({events.mark});
    nm = length(markers);
    
    for i = 1:nm
      
      % Create stimulus group
      % /nirs{i}/stim{j}
      
      nirs_stim_group = create_group(nirs_group, ['stim' num2str(i)]);
      snirf.nirs(gidx).stim(i).name = write_var_string(nirs_stim_group, 'name', markers{i});
      
      marker_idx = find(marker_perm == i);
      ne = length(marker_idx);
      
      stim_data = zeros(ne,3);
      stim_data(:,1) = [events(marker_idx).timestamp];
      stim_data(:,2) = 0;
      stim_data(:,3) = 1;
      snirf.nirs(gidx).stim(i).data = write_double_fixed(nirs_stim_group, 'data', stim_data);
      
      H5G.close(nirs_stim_group);
      
    end
    
  end

  % Create aux groups
  %
  auxi = 1;
  
  if isfield(data(gidx), 'chn_sat')
      % Per frame saturation
      if size(data(gidx).chn_sat, 2) > 1
          nirs_aux_sat = create_group(nirs_group, ['aux' num2str(auxi)]); 
          snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_sat, 'name', 'saturationFlags');
          snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_sat, 'time', [0 data(gidx).chn_dt]);
          snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single_compressed(nirs_aux_sat, 'dataTimeSeries', data(gidx).chn_sat.');
          H5G.close(nirs_aux_sat);
          auxi = auxi+1;
      end
  end
    
  % Temperature
  if isfield(data(gidx), 'node_temp')
      nirs_aux_temp = create_group(nirs_group, ['aux' num2str(auxi)]); 
      snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_temp, 'name', 'temperature');
      snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_temp, 'time', [0 data(gidx).chn_dt]);
      snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single(nirs_aux_temp, 'dataTimeSeries',  data(gidx).node_temp.');
      H5G.close(nirs_aux_temp);
      auxi = auxi+1;    
  end
  
  
  % MPU data
  if isfield(data(gidx), 'node_acc') && isfield(data(gidx), 'node_gyr')
    
    nirs_aux_mpu = create_group(nirs_group, ['aux' num2str(auxi)]); 
    snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_mpu, 'name', 'accel_x');
    snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_mpu, 'time', [0 data(gidx).node_mpu_dt]);
    snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single(nirs_aux_mpu, 'dataTimeSeries',  squeeze(data(gidx).node_acc(:,1,:)).');
    H5G.close(nirs_aux_mpu);
    auxi = auxi+1;    
    
    nirs_aux_mpu = create_group(nirs_group, ['aux' num2str(auxi)]); 
    snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_mpu, 'name', 'accel_y');
    snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_mpu, 'time', [0 data(gidx).node_mpu_dt]);
    snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single(nirs_aux_mpu, 'dataTimeSeries',  squeeze(data(gidx).node_acc(:,2,:)).');
    H5G.close(nirs_aux_mpu);
    auxi = auxi+1;
    
    nirs_aux_mpu = create_group(nirs_group, ['aux' num2str(auxi)]); 
    snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_mpu, 'name', 'accel_z');
    snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_mpu, 'time', [0 data(gidx).node_mpu_dt]);
    snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single(nirs_aux_mpu, 'dataTimeSeries',  squeeze(data(gidx).node_acc(:,3,:)).');
    H5G.close(nirs_aux_mpu);
    auxi = auxi+1;

    nirs_aux_mpu = create_group(nirs_group, ['aux' num2str(auxi)]); 
    snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_mpu, 'name', 'gyro_x');
    snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_mpu, 'time', [0 data(gidx).node_mpu_dt]);
    snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single(nirs_aux_mpu, 'dataTimeSeries',  squeeze(data(gidx).node_gyr(:,1,:)).');
    H5G.close(nirs_aux_mpu);
    auxi = auxi+1;    
    
    nirs_aux_mpu = create_group(nirs_group, ['aux' num2str(auxi)]); 
    snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_mpu, 'name', 'gyro_y');
    snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_mpu, 'time', [0 data(gidx).node_mpu_dt]);
    snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single(nirs_aux_mpu, 'dataTimeSeries',  squeeze(data(gidx).node_gyr(:,2,:)).');
    H5G.close(nirs_aux_mpu);
    auxi = auxi+1;
    
    nirs_aux_mpu = create_group(nirs_group, ['aux' num2str(auxi)]); 
    snirf.nirs(gidx).aux(auxi).name = write_var_string(nirs_aux_mpu, 'name', 'gyro_z');
    snirf.nirs(gidx).aux(auxi).time = write_double(nirs_aux_mpu, 'time', [0 data(gidx).node_mpu_dt]);
    snirf.nirs(gidx).aux(auxi).dataTimeSeries = write_single(nirs_aux_mpu, 'dataTimeSeries',  squeeze(data(gidx).node_gyr(:,3,:)).');
    H5G.close(nirs_aux_mpu);
    auxi = auxi+1;
    
  end
     
  % Close root group
  H5G.close(nirs_group);
    
end
  
  % Close the file, we're done!
  H5F.close(fid);
    
end


function [ml] = write_measlist(nirs_data_group, chn_perm, gi, enum, glch)

  nch = size(glch,1);
  
  % Rather than special case, we'll just write out a standard permutation vector, because
  % we're going to operate over this in a loop anyway.
  if isempty(chn_perm)
    chn_perm = 1:nch;
  end
 
  for ci = 1:nch
    
    nirs_mli_group = create_group(nirs_data_group, ['measurementList' num2str(ci)]);
    
    % glch(ci, 1) -> global source index of channel ci
    % glch(ci, 2) -> global wavelength index of channel ci
    % glch(ci, 3) -> global detector index of channel ci

    ml(ci).sourceIndex = write_int32(nirs_mli_group, 'sourceIndex', glch(chn_perm(ci),1));
    ml(ci).detectorIndex = write_int32(nirs_mli_group, 'detectorIndex', glch(chn_perm(ci),3));
    ml(ci).wavelengthIndex = write_int32(nirs_mli_group, 'wavelengthIndex', glch(chn_perm(ci),2));
    ml(ci).dataType = write_int32(nirs_mli_group, 'dataType', int32(1));
    ml(ci).dataTypeIndex = write_int32(nirs_mli_group, 'dataTypeIndex', int32(1));
    
    % Add source power (this has to be acquired via the canonical enumeration)
    ch = enum.groups(gi).channels(chn_perm(ci));
    src_node_idx = ch.src_node_idx;
    src_node = enum.groups(gi).nodes(src_node_idx);
    src_pwr = src_node.srcs(ch.src_idx).power;
    ml(ci).sourcePower = write_double(nirs_mli_group, 'sourcePower', double(src_pwr));
    
    % Note: it is may be inappropriate to place this information in these fields as we have
    % chosen to export in a global format. We will instead output this in our global to
    % local metadata, which is more consistent with the intent of the specification.
    %
    % write_int32(nirs_mli_group, 'sourceModuleIndex', enum.groups(gi).channels(chn_perm(ci)).src_node_idx);  
    % write_int32(nirs_mli_group, 'detectorModuleIndex', enum.groups(gi).channels(chn_perm(ci)).det_node_idx);
           
  end
  
end

  
function [probe] = write_probe(nirs_probe_group, enum, gidx, glsrc, gldet, glwl)
  
  probe.wavelengths = write_double(nirs_probe_group, 'wavelengths', double(glwl));
   
  nwl = length(glwl);
  nsrc = size(glsrc,2);
  ndet = size(gldet,2);
  
  sourcePos2D = zeros(nsrc, 2);
  sourcePos3D = zeros(nsrc, 3);
  sourceLabels = cell(nsrc, nwl);

  % Buidl source positions
  for i = 1:nsrc

    node_idx = glsrc(i).node_idx;
    node_id = enum.groups(gidx).nodes(node_idx).id;
    optode_idx = glsrc(i).optode_idx;
    dock_idx = enum.groups(gidx).layout.dockmap(node_id);       

    node_optode = enum.groups(gidx).nodes(node_idx).optodes(optode_idx);
    optode = enum.groups(gidx).layout.docks(dock_idx).optodes(optode_idx);

    assert(node_optode.name == optode.name);
    
    sourcePos2D(i,:) =  [optode.coords_2d.x optode.coords_2d.y];
    sourcePos3D(i,:) =  [optode.coords_3d.x optode.coords_3d.y optode.coords_3d.z];
    
    for j = 1:nwl
      sourceLabels{i,j} = ['N' num2str(node_id) '-' optode.name num2str(glwl(j))];
    end
    
  end
  
  probe.soucePos2D = write_double(nirs_probe_group, 'sourcePos2D', sourcePos2D);
  probe.soucePos3D = write_double(nirs_probe_group, 'sourcePos3D', sourcePos3D);
  probe.souceLabels = write_var_string(nirs_probe_group, 'sourceLabels', sourceLabels);
  
  % Build detector positions
  detectorPos2D = zeros(ndet, 2);
  detectorPos3D = zeros(ndet, 3);
  detectorLabels = cell(ndet, 1);
  
  for i = 1:ndet

    node_idx = gldet(i).node_idx;
    node_id = enum.groups(gidx).nodes(node_idx).id;
    optode_idx = gldet(i).optode_idx;
    dock_idx = enum.groups(gidx).layout.dockmap(node_id);       

    node_optode = enum.groups(gidx).nodes(node_idx).optodes(optode_idx);
    optode = enum.groups(gidx).layout.docks(dock_idx).optodes(optode_idx);

    assert(node_optode.name == optode.name);
    
    detectorPos2D(i,:) =  [optode.coords_2d.x optode.coords_2d.y];
    detectorPos3D(i,:) =  [optode.coords_3d.x optode.coords_3d.y optode.coords_3d.z];
    detectorLabels{i} = ['N' num2str(node_id) '-' optode.name];
    
  end
      
  probe.detectorPos2D = write_double(nirs_probe_group, 'detectorPos2D', detectorPos2D);
  probe.detectorPos3D = write_double(nirs_probe_group, 'detectorPos3D', detectorPos3D);
  probe.detectorLabels = write_var_string(nirs_probe_group, 'detectorLabels', detectorLabels);
  
  % Write landmarks  
  if ~isempty(enum.groups(gidx).layout.landmarks)
    
    nl = length(enum.groups(gidx).layout.landmarks);
    landmarkLabels = cell(nl, 1);
    landmarkPos3D = zeros(nl, 4);
 
    for i = 1:nl
      landmark = enum.groups(gidx).layout.landmarks(i);
      landmarkLabels{i} = landmark.name;
      landmarkPos3D(i, 1:3) = [landmark.coords_3d.x landmark.coords_3d.y landmark.coords_3d.z];
      landmarkPos3D(i, 4) = i;
    end
        
    probe.landmarkPos3D = write_double(nirs_probe_group, 'landmarkPos3D', double(landmarkPos3D));
    probe.landmarkLabels = write_var_string(nirs_probe_group, 'landmarkLabels', landmarkLabels);
    
  end
 
end

function [gid] = create_group(base, path)
pl = 'H5P_DEFAULT';
gid = H5G.create(base, path, pl, pl, pl);
end

%%% HDF5 writing functions
%
% Notes:
%
% - When a single element is provided, a scalar dataset is written
% - When a vector (either row or column) is provided, a rank 1 dataset is written
% - When a matrix is provided, it is transposed before write

function data = write_chn_dat_block(nirs_data_group, chn_perm, data)

h5_data_size = fliplr(size(data));
h5_data_rank = ndims(data);

tp = H5T.copy('H5T_IEEE_F32LE');                                  % Type desc.
ds = H5S.create_simple(h5_data_rank, h5_data_size, h5_data_size); % Dataspace
pl = H5P.create('H5P_DATASET_CREATE');                            % Porperty list

% For now (e.g. before lazy loading and channel by channel manipulation), we will write the
% data time series without chunking and compression (typical ratios appear to be around 1.2
% when chunking along the time axis).
%
% h5_chunk_size = fliplr([size(data,1) 1]);
% h5_chunk_size = fliplr([1 size(data,2)]);
% h5_dflt_lvl = 7;
%
% H5P.set_chunk(pl, h5_chunk_size);           % Set chunk size to be one frame
% H5P.set_deflate(pl, h5_dflt_lvl);           % Set deflate compression

% Dataset
dset = H5D.create(nirs_data_group, 'dataTimeSeries', tp, ds, pl);

% Write
if ~isempty(chn_perm)
  H5D.write(dset, tp, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data(chn_perm,:));
else
  H5D.write(dset, tp, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
end

% Clean up
H5P.close(pl);
H5T.close(tp);
H5S.close(ds);
H5D.close(dset);

end



function str = write_var_string(base, path, str)
% Write a string or a cell array of strings as a scalar, vector or matrix of variable length
% HDF5 strings.

% Prepare
tp = H5T.copy('H5T_C_S1');                  % String type
H5T.set_size(tp, 'H5T_VARIABLE');           % Set string variable length

if ~iscell(str)
  
  ds = H5S.create('H5S_SCALAR');   
  pl = H5P.create('H5P_DATASET_CREATE');      % Porperty list
  dset = H5D.create(base, path, tp, ds, pl);  % Dataset

  % Write
  H5D.write(dset, tp, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', str);
  
else
  
  if rankn(str) == 1
    % Vector
    len = numel(str);
    ds = H5S.create_simple(1, len, len);
  elseif rankn(str) == 2
    % Matrix (must transpose)
    str = str.';
    ds = H5S.create_simple(2, fliplr(size(str)), fliplr(size(str)));
  else
    error('Data dimensionality not supported');
  end
   
  pl = H5P.create('H5P_DATASET_CREATE');      % Porperty list
  dset = H5D.create(base, path, tp, ds, pl);  % Dataset

  % Write 
  H5D.write(dset, tp, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', str);
    
  
end

% Clean up
H5P.close(pl);
H5T.close(tp);
H5S.close(ds);
H5D.close(dset);

end

function val = write_int32(base, path, val)
tp = H5T.copy('H5T_STD_I32LE');
write_core(base, path, int32(val), tp, false, false);
end

function val = write_int32_compressed(base, path, val)
tp = H5T.copy('H5T_STD_I32LE');
write_core(base, path, int32(val), tp, false, true);
end

function val = write_int8(base, path, val)
tp = H5T.copy('H5T_STD_I8LE');
write_core(base, path, int8(val), tp, false, false);
end

function val = write_double(base, path, val)
tp = H5T.copy('H5T_IEEE_F64LE');
write_core(base, path, double(val), tp, false, false);
end

function val = write_double_fixed(base, path, val)
tp = H5T.copy('H5T_IEEE_F64LE');
write_core(base, path, double(val), tp, true, false);
end

function val = write_single(base, path, val)
tp = H5T.copy('H5T_IEEE_F32LE');
write_core(base, path, single(val), tp, false, false);
end

function val = write_single_compressed(base, path, val)
tp = H5T.copy('H5T_IEEE_F32LE');
write_core(base, path, single(val), tp, false, true);
end

function write_core(base, path, val, tp, fix_shape, compress)

% Dataspace
if fix_shape
  val = val.';
  ds = H5S.create_simple(2, fliplr(size(val)), fliplr(size(val)));
  h5_chunk_size = fliplr(size(val));
  h5_chunk_size(1:(end-1)) = 1;
else
  if isscalar(val)
    % Scalar value
    ds = H5S.create('H5S_SCALAR'); 
    h5_chunk_size = 1;
    elseif rankn(val) == 1
    % Vector
    len = numel(val);
    ds = H5S.create_simple(1, len, len);
    h5_chunk_size = len;
  elseif rankn(val) == 2
    % Matrix (must transpose)
    val = val.';
    ds = H5S.create_simple(2, fliplr(size(val)), fliplr(size(val)));
    h5_chunk_size = fliplr(size(val));
    h5_chunk_size(1:(end-1)) = 1;
  else
    H5T.close(tp);
    error('Data dimensionality not supported');
  end
end

pl = H5P.create('H5P_DATASET_CREATE');      % Porperty list

if compress
  
  % If compression is required, we will chunk by the fastest varying dimension in the
  % dataset (the last dimension in HDF5).
  %
  h5_dflt_lvl = 7;
  H5P.set_chunk(pl, h5_chunk_size);           % Set chunk size to be one frame
  H5P.set_deflate(pl, h5_dflt_lvl);           % Set deflate compression
  
  fprintf('Applying compression with chunk size: [');
  fprintf('%g ', h5_chunk_size);
  fprintf(']\n');

end

dset = H5D.create(base, path, tp, ds, pl);  % Dataset

% Write
H5D.write(dset, tp, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', val);

% Clean up
H5P.close(pl);
H5S.close(ds);
H5D.close(dset);
H5T.close(tp);

end

function r = rankn(d)
  if numel(d) == 1
    r = 1;
  else
    r = sum(size(d) > 1);
  end  
end



