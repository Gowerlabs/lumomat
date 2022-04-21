function validate_layout(enum)
%VALIDATE_LAYOUT Check consistency of the enumeration and the layout

% We only handle one group for now.
assert(length(enum.groups) == 1);

idl = hex2dec(enum.groups.layout.id);
idg = hex2dec(enum.groups.id);

% Check that the reocrded IDs match up
if idl ~= idg
  warning('Layout ID (%s) does not match group ID (%s), enumeration may be inconsistent',...
          enum.groups.layout.id, enum.groups.id);     
end

% Check that all node IDs have docks
nids = [enum.groups.nodes.id];
dids = [enum.groups.layout.docks.id];

if ~all(ismember(nids, dids))
  error('Layout ID (%s) does not contain dock positions for every node in the enumeration', ...
        enum.groups.layout.id);
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

