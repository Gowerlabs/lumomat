function write_SNIRF(snirffn,  enum, data, events, varargin)
% LUMOFILE.WRITE_SNIRF Write LUMO data to SNIRF format
%
%   LUMOFILE.WRITE_SNIRF(snirffn, enum, data, events)
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
%   'draft_chdesc':         Use the draft modifications to the channel descriptor groups
%                           which significant reduce file size and improve performance, as
%                           discussed here: https://github.com/fNIRS/snirf/issues/103.
%
%   Details:
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
%   The following LUMO specific information is written to the file:
%
%
% See also LUMO_READ
%
%
%   (C) Gowerlabs Ltd., 2022
%

%%% TODO
%
% Add details of the additional LUMO specific fields
% Add global saturation flags
% Add temporal flags
%

% Check that a layout is available for every group
ng = length(enum.groups);

for gidx = 1:ng
  if isempty(enum.groups(gidx).layout)
    error('Input LUMO enumeration does not contain an embedded layout file');
  end
end

% try
%     if(isa(fname,'H5ML.id'))
%         fid=fname;
%     else
%         fid = H5F.create(fname, 'H5F_ACC_TRUNC', H5P.create('H5P_FILE_CREATE'), H5P.create('H5P_FILE_ACCESS'));
%     end
%     obj2h5(rootname,data,fid,1,opt);
% catch ME
%     if(exist('fid','var') && fid>0)
%         H5F.close(fid);
%     end
%     rethrow(ME);
% end

fprintf('Writing SNIRF file %s...\n', snirffn);

% Open the file
try
  fid = H5F.create(snirffn, 'H5F_ACC_TRUNC', H5P.create('H5P_FILE_CREATE'), H5P.create('H5P_FILE_ACCESS'));
catch
  fprintf('Error creating SNIRF file %s\n', snirffn);
  rethrow(e);
end

% Write format version
%
write_var_string(fid, '/formatVersion', '1.0');

% Over each group
%
for gidx = 1:ng
  
  % We use a single data block
  bi = 1;
  
  % Build the global spectroscopic mapping
  %
  try
    [glch, glsrc, gldet, glwl] = lumofile.map_gs(enum, 'group', gidx);
  catch e
    fprintf('Error forming gloabl spectroscopic mapping\n');
    rethrow(e)
  end
  
  % Create NIRS root
  % /nirs{i}
  %
  nirs_group = create_group(fid, ['nirs' num2str(gidx)]);
  
  % Create metadata
  % /nirs{i}/metaDataTags
  %
  
  % Required fields
  nirs_meta_group = create_group(nirs_group, 'metaDataTags');
  write_var_string(nirs_meta_group, 'SubjectID', 'Subject Unknown');
  write_var_string(nirs_meta_group, 'MeasurementDate', 'unknown');
  write_var_string(nirs_meta_group, 'MeasurementTime', 'unknown');
  write_var_string(nirs_meta_group, 'LengthUnit', 'mm');
  write_var_string(nirs_meta_group, 'TimeUnit', 'ms');
  write_var_string(nirs_meta_group, 'FrequencyUnit', 'Hz');
  
  % LUMO sepcific fields
  write_var_string(nirs_meta_group, 'ManufacturerName', 'Gowerlabs');
  write_var_string(nirs_meta_group, 'Model', 'LUMO');

  % write_var_string(nirs_meta_group, 'Hub SN');
  % write_var_string(nirs_meta_group, 'Layout', jsonencode...)
  % write_var_string(nirs_meta_group, 'Nodes');
  % write_var_string(nirs_meta_group, 'Group UID')
   
  
  H5G.close(nirs_meta_group);
  
  % Create data block
  % /nirs{i}/data{i}
  %
  nirs_data_group = create_group(nirs_group, ['data' num2str(bi)]);
  write_chn_dat_block(nirs_data_group, data(gidx).chn_dat);       % Write data /nirs{i}/data1
  write_double(nirs_data_group, 'time', [0 data(gidx).chn_dt]);   % Write /nirs{i}/time
  write_measlist(nirs_data_group, gidx, enum, glch);              % Write measurementList 
  H5G.close(nirs_data_group);
  
  % Create probe block
  % /nirs{i}/probe
  %
  nirs_probe_group = create_group(nirs_group, 'probe');
  write_probe(nirs_probe_group, enum, gidx, glsrc, gldet, glwl);
  H5G.close(nirs_probe_group);
  
  % Create stimulus group
  %
  warning('Not writing stimulus');

  % Create aux group
  %
  %%% TODO
  %
  % Add MPU
  % Add temperature
  % Add saturation
  warning('Not writing aux data');
  

    
  % Close root group
  H5G.close(nirs_group);
    
end
  
  % Close the file, we're done!
  H5F.close(fid);
    
end


function write_chn_dat_block(nirs_data_group, data)

h5_data_size = fliplr(size(data));
h5_data_rank = ndims(data);
h5_chunk_size = fliplr([size(data,1) 1]);
h5_dflt_lvl = 3;

tp = H5T.copy('H5T_IEEE_F32LE');            % Single precision little endien
% Dataspace
ds = H5S.create_simple(h5_data_rank, h5_data_size, h5_data_size);
pl = H5P.create('H5P_DATASET_CREATE');      % Porperty list
H5P.set_chunk(pl, h5_chunk_size);           % Set chunk size to be one frame
H5P.set_deflate(pl, h5_dflt_lvl);           % Set deflate compression
% Dataset
dset = H5D.create(nirs_data_group, 'dataTimeSeries', tp, ds, pl);

% Write
H5D.write(dset, tp, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);

% Clean up
H5P.close(pl);
H5T.close(tp);
H5S.close(ds);
H5D.close(dset);

warning("NOT WRITING TIME");

end

function write_measlist(nirs_data_group, gi, enum, glch)

  nch = size(glch,1);
 
  for ci = 1:nch
    
    nirs_mli_group = create_group(nirs_data_group, ['measurementList' num2str(ci)]);
    
    % glch(ci, 1) -> global source index of channel ci
    % glch(ci, 2) -> global wavelength index of channel ci
    % glch(ci, 3) -> global detector index of channel ci

    write_int32(nirs_mli_group, 'sourceIndex', glch(ci,1))
    write_int32(nirs_mli_group, 'detectorIndex', glch(ci,3))
    write_int32(nirs_mli_group, 'wavelengthIndex', glch(ci,2))
    write_int32(nirs_mli_group, 'dataType', int32(1))
    write_int32(nirs_mli_group, 'dataTypeIndex', int32(1))
    
    %% TODO can we include this field even if we use global indexing? 
    %
    write_int32(nirs_mli_group, 'sourceModuleIndex', enum.groups(gi).channels(ci).src_node_idx);  
    write_int32(nirs_mli_group, 'detectorModuleIndex', enum.groups(gi).channels(ci).det_node_idx);
    
    if ci == 1
      warning('Not writing source power');
    end
    %write_float(nirs_ml_group, 'srcPower', X)
       
  end
  
end

  
function write_probe(nirs_probe_group, enum, gidx, glsrc, gldet, glwl)
  
  write_double(nirs_probe_group, 'wavelengths', glwl);
   
  nwl = length(glwl);
  nsrc = size(glsrc,2);
  ndet = size(gldet,2);
  
  sourcePos2D = zeros(nsrc, 2);
  sourcePos3D = zeros(nsrc, 3);
  sourceLabels = cell(nsrc, nwl);

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
  
  write_double(nirs_probe_group, 'sourcePos2D', sourcePos2D);
  write_double(nirs_probe_group, 'sourcePos3D', sourcePos3D);
  
    % Build source positions
    
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
      
  write_double(nirs_probe_group, 'detectorPos2D', detectorPos2D);
  write_double(nirs_probe_group, 'detectorPos3D', detectorPos3D);

  %% TODO
  %
  % In the global scehme we can add some additional data for the user here. If we are
  % permitted to store the module indices despite using a global enumeration, these can
  % simply be 'A', '!', etc.
  %
  % sourceLabels
  % detectorLabels
  
  %% TODO
  %
  % Grab data from the layout
  %
  %write_double(nirs_prove_group, 'landmarkPos3D', X);
  %write_double(nirs_probe_group, 'landmarkLabels', X);
  
 
end

function [gid] = create_group(base, path)
pl = 'H5P_DEFAULT';
gid = H5G.create(base, path, pl, pl, pl);
end

function write_var_string(base, path, str)

% Prepare
tp = H5T.copy('H5T_C_S1');                  % String type
H5T.set_size(tp, 'H5T_VARIABLE');           % Set string variable length
ds = H5S.create('H5S_SCALAR');   
pl = H5P.create('H5P_DATASET_CREATE');      % Porperty list
dset = H5D.create(base, path, tp, ds, pl);  % Dataset

% Write
H5D.write(dset, tp, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', str);

% Clean up
H5P.close(pl);
H5T.close(tp);
H5S.close(ds);
H5D.close(dset);

end

function write_int32(base, path, val)
tp = H5T.copy('H5T_STD_I32LE');
write_core(base, path, int32(val), tp);
end

function write_double(base, path, val)
tp = H5T.copy('H5T_IEEE_F64BE');
write_core(base, path, double(val), tp);
end

function write_single(base, path, val)
tp = H5T.copy('H5T_IEEE_F32BE');
write_core(base, path, single(val), tp);
end

function write_core(base, path, val, tp)

% Dataspace
if isscalar(val)
  % Scalar value
  ds = H5S.create('H5S_SCALAR');     
elseif rankn(val) == 1
  % Vector
  len = numel(val);
  ds = H5S.create_simple(1, len, len);
elseif rankn(val) == 2
  % Matrix
  ds = H5S.create_simple(2, fliplr(size(val)), fliplr(size(val)));
else
  H5T.close(tp);
  error('Data dimensionality not supported');
end
    
pl = H5P.create('H5P_DATASET_CREATE');      % Porperty list
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



