function [layout] = proc_layout(layout_raw)
% LUMOFILE.PROC_LAYOUT Process a layout in JSON format
%
% [layout] = LUMOFILE.PROC_LAYOUT(layout_raw) 
%
% LUMOFILE.PROC_LAYOUT is a helper function typically called by LUMOFILE.READ_LAYOUT. 
%
%   (C) Gowerlabs Ltd., 2022
%

% Get the variables
try
  lf_lo_group_uid = layout_raw.group_uid;
  lf_lo_dims_2d = layout_raw.dimensions.dimensions_2d;
  lf_lo_dims_3d = layout_raw.dimensions.dimensions_3d;
  
  lf_lo_landmarks = optfieldci(layout_raw, 'landmarks'); % Case varies, accept either
  
  if ~isempty(lf_lo_landmarks)
    assert(isfield(lf_lo_landmarks, 'name'));
    assert(isfield(lf_lo_landmarks, 'x'));
    assert(isfield(lf_lo_landmarks, 'y'));
    assert(isfield(lf_lo_landmarks, 'z'));
    
    % Restructure the layout fields to permit 2D/3D representations
    for li = 1:length(lf_lo_landmarks)
      lf_lo_landmarks(li).coords_3d.x = lf_lo_landmarks(li).x;
      lf_lo_landmarks(li).coords_3d.y = lf_lo_landmarks(li).y;
      lf_lo_landmarks(li).coords_3d.z = lf_lo_landmarks(li).z;
    end
    lf_lo_landmarks = rmfield(lf_lo_landmarks, {'x', 'y', 'z'});
    
  end    

  % Over each dock (first pass for sorting)
  nd = length(layout_raw.docks);
  dockidsort = zeros(nd, 1);
  for di = 1:nd
    docknum = strsplit(layout_raw.docks(di).dock_id, '_');
    dockidsort(di) = str2num(docknum{2});
  end
  [~, dock_perm] = sort(dockidsort);
  
  % Over each node(second pass for construction)
  for dii = 1:nd
    
    di = dock_perm(dii);
    
    assert(length(layout_raw.docks(di).optodes) == 7);
    for oi = 1:7
      [cs_optode, opt_idx] = trans_optode_desc(layout_raw.docks(di).optodes(oi));
      optperm(oi) = opt_idx;
      optodes(oi) = cs_optode;
    end
    
    optodes = optodes(optperm);
    
    lf_lo_docks(dii) = struct('id', dii, 'optodes', optodes);
  end
  
catch e
  fprintf('LUMO file (%s): error parsing layout structure from file\n');
  rethrow(e);
end

% Form the ID to index map for later convenience
dockids = [lf_lo_docks.id];
maxdid = max(dockids);
dockmap = zeros(maxdid, 1);
dockmap(dockids) = 1:length(dockids);

[uid_hex, uid_name] = lumomat.norm_gid(lf_lo_group_uid);

layout = struct('id', uid_hex, ...
  'name', uid_name, ...
  'dims_2d', lf_lo_dims_2d, ...
  'dims_3d', lf_lo_dims_3d, ...
  'landmarks', lf_lo_landmarks, ...
  'docks', lf_lo_docks,...
  'dockmap', dockmap);

% lumofile.validate_layout(layout);

end


% trans_optode_desc
%
% Translate the .lumo template layout format into the canonical format.
function [cs_optode, opt_idx] = trans_optode_desc(lf_optode)

switch lf_optode.optode_id
  case 'optode_1'
    opt_idx = 1;
    opt_name = '1';
  case 'optode_2'
    opt_idx = 2;
    opt_name = '2';
  case 'optode_3'
    opt_idx = 3;
    opt_name = '3';
  case 'optode_4'
    opt_idx = 4;
    opt_name = '4';
  case 'optode_a'
    opt_idx = 5;
    opt_name = 'A';
  case 'optode_b'
    opt_idx = 6;
    opt_name = 'B';
  case 'optode_c'
    opt_idx = 7;
    opt_name = 'C';
  otherwise
    error('Error parsing optode structure (optode id %s)', lf_optode.optode_id);
end

coords_2d = lf_optode.coordinates_2d;
coords_3d = lf_optode.coordinates_3d;

cs_optode = struct(...
  'name', opt_name, ...
  'coords_2d', coords_2d,...
  'coords_3d', coords_3d);

end
