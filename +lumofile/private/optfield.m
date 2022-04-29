function fieldval = optfield(s, name)
% OPTFIELD
%
% Attempt to get an optional field from a structure, returning an empty array if unavailable.
%

if isfield(s, name)
  fieldval = getfield(s, name);
else
  fieldval = [];
end

end