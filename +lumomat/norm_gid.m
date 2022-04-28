function [uid_hex, uid_name] = norm_gid(uid)
%NORM_GID Normalise an input group UID to hex string and name
%
% Input must be numeric or a corresponding hexadecimal string

   if(uid > (flintmax-1))
        warning('Group ID out of range, recording may be invalid');
   end

  % Normalise the UID to a hex string
  if isnumeric(uid)
    uid =  ['0x' lower(dec2hex(uid, 8))];
  end
  uid_hex = uid;

  % Now form the name by converting UID to string
  uid = str2num(uid);

  if(uid > 20155392)
    if(uid > (flintmax-1))
      group_name = '';
    else
      group_name = dec2base(uid, 36);
    end
  else
    group_name = sprintf('GA%05d', uid);
  end 
  uid_name = group_name;

end
