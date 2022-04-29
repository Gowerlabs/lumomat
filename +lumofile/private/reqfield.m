% reqfield
%
% Get a required field from a structure, throwing an appropriate error if unavailable.
%
function fieldval = reqfield(s, name, lumodir_fn)

if isfield(s, name)
  fieldval = getfield(s, name);
else
  error('LUMO file (%s): required field %s missing from structure', lumodir_fn, name)
end

end
