function [layout_file_dir, layout_file_name] = check_app_data_for_layout(group_id)

%% Get Directory
local_app_data_environment = getenv('localappdata');

% Check if localappdata directory exists
if(isempty(local_app_data_environment))
    error("Unable to retrieve layout files: This Operating System does not have a localappdata environment.");
end

% Get Directory with layout files
layout_file_dir = fullfile(local_app_data_environment, "Gowerlabs", "Lumo");

% Check if directory exists
if(not(isfolder(layout_file_dir)))
    error("Unable to retrieve layout files: This machine either doesn't have the LumoView software installed, or was never had it used before.");
end

