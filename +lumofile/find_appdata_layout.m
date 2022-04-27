function [layout_file_path] = find_appdata_layout(group_id)
% LUMOFILE.find_appdata_layout finds the apropriate appData if it exists
% and returns it's full path as a string, otherwise it returns an empty
% matrix.

% Group_id refers to the id of the cap, this value can be found on the cap
% and used as an input for this function.

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

%% Get group_id

if(isnumeric(group_id))
    % Simple if value is numeric
    group_id_num = group_id;
else
    %Otherwise, checking if value is a cell of char or char.
    if(isa(group_id,'char'))
        group_id_str = group_id;
    elseif(isa(group_id,'string'))
        group_id_str = group_id;
    elseif(isa(group_id,'cell'))
        if(isa(group_id{1}, 'char'))
            group_id_str = group_id{1};
        elseif(isa(group_id{1},'string'))
            group_id_str = group_id;
        else
            error("Unable to retrieve layout files: group_id is not a char, cell with a char or numeric type.");
        end
    else
        error("Unable to retrieve layout files: group_id is not a char, cell with a char or numeric type.");
    end
    
    %Checking if char contains hex or dec value.
    
    try
        first2Chars = extractBetween(group_id_str, 1, 2) == "0x"
        
        if(first2Chars == "0x")
        %% Group ID is in hexadecimal
            group_id_str_arr = split(group_id_str,"x");
            group_id_num = hex2dec(group_id_str_arr{2});
             
        elseif(first2Chars == "GA")
        %% Group ID is using V2 Cap IDs
            group_id_str_arr = split(group_id_str,"A");
            group_id_num = group_id_str_arr{2};
             
        elseif(first2Chars == "C0")
        %% Group ID is using V1 cap IDs
            group_id_num = base2dec{group_id_str};
            
        else
            %Check if valid number
            str2num(group_id_str);
            group_id_num = group_id_str;
        end
    catch
        error("Unable to retrieve layout files: group_id contains invalid characters.");
    end
    
    %Check for precision, because this value can be any UInt64 value, but
    %double precision floats loose their precision earlier and num2hex
    %uses doubles this may cause issues and the user should be warned.
    if(group_id_num == group_id_num +1)
        warning("group_id's value is too large to be accuratley stored as a " + class(group_id) + ", the produced layout file may be incorrect");
    end
end
    
%% Get File name

layout_file_name = "coordinates_" + group_id_num + ".json";

if(not(isfile(fullfile(layout_file_dir, layout_file_name))))
    error("Unable to retrieve layout files: file does not exist on this computer.");
end

layout_file_path = fullfile(layout_file_dir, layout_file_name);
