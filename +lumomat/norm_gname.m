function [group_id] = norm_gname(name)
%NORM_GNAME Normalise an input group name to a numeric group ID
%
% Input must be a group name such as 'CA0001', or 'GA0001'

% Normalise strings
name = convertStringsToChars(name);

if ~ischar(name)
  error('Input group name must be a string');
end

if length(name) < 5
  error('Group name too short');
end

name_prefix = name(1:2);

if strcmp(name_prefix, 'GA')
  
    % Group ID is using V2 Cap IDs
    group_id = str2num(name(3:end));
        
elseif (strcmp(name_prefix, 'C0') || strcmp(name_prefix, 'CA'))
    
    % Group ID is using V1 cap IDs
    group_id = base2dec(name,36);
  
else
    error('Group name unknown')
end
    
end


