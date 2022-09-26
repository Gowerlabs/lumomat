function validate_layout(enum, gidx)
%VALIDATE_LAYOUT Check consistency of the enumeration and the layout

idl = sscanf(enum.groups(gidx).layout.id, '%x');
idg = sscanf(enum.groups(gidx).id, '%x');

% Check that the reocrded IDs match up
if idl ~= idg
  warning('Layout ID (%s) does not match group ID (%s), enumeration may be inconsistent',...
          enum.groups(gidx).layout.id, enum.groups(gidx).id);     
end

% Check that all node IDs have docks
nids = [enum.groups(gidx).nodes.id];
dids = [enum.groups(gidx).layout.docks.id];

if ~all(ismember(nids, dids))
  error('Layout ID (%s) does not contain dock positions for every node in the enumeration', ...
        enum.groups(gidx).layout.id);
end

% % Check every optode can be accessed
% for i = 1:length(nids)
%   
%   
%   
%   
% end
%   

end

