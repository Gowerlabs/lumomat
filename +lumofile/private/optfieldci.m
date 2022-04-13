function [val] = optfieldci(s, name)
%OPTFIELDCI Return optional case insensitive field from structure
%
% [val] = optfieldci(s, name)
%
%
names   = fieldnames(s);
isField = strcmpi(name,names);

if any(isField)
  val = s.(names{isField});
else
  val = [];
end

end


