function [layout_file_path] = find_install_layout(group_id)
% LUMOFILE.FIND_INSTALL_LAYOUT Locate layout file from LUMOview installation
%
%   [layout_file_path] = FIND_INSTALL_LAYOUT(group_id)
%
% FIND_INSTALL_LAYOUT attempts to locate a layout file for the specified group ID in the
% standard installation paths used by the LUMOview software.
%
%   Parameters
%
%   group_id:           numeric group ID
%
%   Returns:
%
%   layout_file_path:   a string containing the path to the appropriate file, or an empty
%                       matrix if no suitable file is located.
%
%
%   (C) Gowerlabs Ltd., 2022
%

layout_file_path = [];

    
% Search for the installation layout path
appdata_env = getenv('localappdata');
if isempty(appdata_env)
    return;
end

layout_dir = fullfile(appdata_env, 'Gowerlabs', 'Lumo');

if ~isfolder(layout_dir)
    return;
end

% Get group_id
[uid_hex, uid_name] = lumomat.norm_gid(group_id);
group_id_num = hex2dec(uid_hex);

% Check for file
layout_fn = ['coordinates_' num2str(group_id_num) '.json'];
layout_path_check = fullfile(layout_dir, layout_fn);

if isfile(layout_path_check)
  layout_file_path = layout_path_check;
end


end