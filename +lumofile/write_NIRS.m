function [nirs] = write_NIRS(nirsfn, enum, data, events, varargin)
% LUMOFILE.WRITE_NIRS Covert LUMO data to NIRS format and write to disk
%
%   [nirs] = LUMOFILE.WRITE_NIRS(nirsfn, enum, data, events)
%
%   LUMOFILE.WRITE_NIRS constructs and writes to disk a NIRS format data structure as used 
%   by the Homer2 analysis software. The format of the NIRS structure is detailed in the
%   Homer2 user guide:
%
%   https://www.nmr.mgh.harvard.edu/martinos/software/homer/HOMER2_UsersGuide_121129.pdf
%
%   Note: there is some ambiguity as to the form of the SD structure generated as part of
%   the NIRS format. This function provides different options for this structure.
%
%   Paramters:
%
%     nirsfn:               The file name of the output NIRS file, if empty, the dta is only
%                           returned and not written to disk.
%
%     enum, data, events:   Data structures returned by LUMOFILE.READ
%
%   Optional Parameters:
%
%   'sd_style': A string describing the style of the SD field to be built:
%
%               'standard':  (default) construct SD output as per NIRS specification with 3D
%                            locations for the sources and detectors.
%
%               'flat':      construct the SD output with a flattened 2D layout of the
%                            source/detector layout. Include a second structure called SD3D
%                            which contains the full three dimensional layout. This format
%                            is used by the DOT-HUB toolbox.
%
%   'group':    An integer specifiying the group index of the enumeration to map. Defaults to
%               group index 1.
%   Returns:
%
%     nirs:     The NIRS format structure.
%
%   Details:
%
%   To construct a NIRS format description of the system, the LUMO enumeration is
%   transformed from the canonical node local format to a globally indexed spectroscopic
%   format, such that channels are described by the tuple:
%
%   (src_position(i), det_postion(j), src_wavelength(k))
%
%   The canonical LUMO format permits more complex representations of the system than can be
%   accommodated by this format, so the actual enumeration is checked to ensure that it can
%   conform to this representation.
%
%   When constructing the global enumeration, the template layout (e.g. collection of all
%   docks in a group/cap) is collapsed to form a description containing only the occupied
%   docks.
%
%   The following additional fields are provided in NIRS structures output from this
%   function, which are not present in the specification:
%
%   SD.SpatialUnit:     A string representing the SI units of distance
%
%   lumo.chn_sat:       Vector (or matrix) of logical values which indicate if a particular
%                       channel is saturated during the recording. When this field is a 
%                       matrix, the saturation flag is provided for every frame in the
%                       recording, for more granular exclusion of data.
%
%   lumo.src_pwr:       Vector of source powers per channel.
%
%   lumo.chn_sort_perm: Permutation which maps from the re-indexed channel listing
%                       constructed for NIRS back to the canonical indexing
%
%
% See also LUMO_READ
%
%
%   (C) Gowerlabs Ltd., 2022
%
%   Portions of code reused from DOT-HUB toolbox with permission.
%

%%% TODOS
%
% - Check with RJC regarding event data construction
%

nirs = lumo_build_NIRS(enum, data, events, varargin{:});

if ~isempty(nirsfn)
  fprintf('Writing NIRS file %s... ', nirsfn);
  save(nirsfn, '-struct', 'nirs', '-v7.3');
  fprintf('complete.\n');
end

end

function [nirs] = lumo_build_NIRS(enum, data, events, varargin)

p = inputParser;
expected_styles = {'standard', 'flat', 'dhtoolbox'};
addOptional(p, 'sd_style', 'standard', @(x) any(validatestring(x, expected_styles)));
addOptional(p, 'group', 1, @(x) (isnumeric(x) && x > 0));
parse(p, varargin{:})

gi = p.Results.group;
sdstyle = p.Results.sd_style;

% Check the group index is sensible
if (gi > length(enum.groups)) || (gi > length(data))
  error('Requested group %d exceeds that available in the enumeration/data',gi)
end

% Check that a layout is available
if isempty(enum.groups(gi).layout)
  error('Input LUMO enumeration does not contain an embedded layout file');
end

% Time vector
t = ((0:(data(gi).nframes - 1))./data(gi).chn_fps).';   % Time [nt x 1]

% Build the global spectroscopic mapping
[glch, glsrc, gldet, glwl] = lumofile.map_gs(enum);

SD.nSrcs = size(glsrc,2);
SD.nDets = size(gldet,2);
SD.Lambda = glwl;

% Build SD mapping
%
% Since LUMO is a modular system it is possible that not all available docks will be
% populated with tiles. The SD mapping is constructed based only on the populated tiles, so
% we loop over the global enumeration and construct the list of optodes accordingly.
%
% There are two styles of mapping available:
%
% 'standard':   SD structure constructed as per specification
%
% 'flat':       SD structure built with two-dimensional flat layout, and an additional SD3D
%               structure is provided with the real 3D layout (required by DOT-HUB toolbox)
%
SD.SrcPos = zeros(SD.nSrcs, 3);
SD.DetPos = zeros(SD.nDets, 3);
SD.SrcPowers = zeros(SD.nSrcs, 2);
SD.SpatialUnit = 'mm';          % Out of spec but required by DHT

if strcmp(sdstyle,'flat')
  SD3D = SD;
end

layout = enum.groups(gi).layout;

% Add source entries
for qi = 1:SD.nSrcs
  
  nidx = glsrc(qi).node_idx;
  oidx = glsrc(qi).optode_idx;
  nid = enum.groups(gi).nodes(nidx).id;
  didx = layout.dockmap(nid);
  
 % Get the nodes sources, find those on the optode
 nd_srcs = enum.groups.nodes(nidx).srcs;
 nd_oidx = [nd_srcs.optode_idx] == oidx;
 assert(all([nd_srcs(nd_oidx).wl] == [735 850]));
 SD.SrcPowers(qi,:) = [nd_srcs(nd_oidx).power];
  
  switch sdstyle
    case 'standard'
      SD.SrcPos(qi, 1) = layout.docks(didx).optodes(oidx).coords_3d.x;
      SD.SrcPos(qi, 2) = layout.docks(didx).optodes(oidx).coords_3d.y;
      SD.SrcPos(qi, 3) = layout.docks(didx).optodes(oidx).coords_3d.z;
      
    case 'flat'
      SD.SrcPos(qi, 1) = layout.docks(didx).optodes(oidx).coords_2d.x;
      SD.SrcPos(qi, 2) = layout.docks(didx).optodes(oidx).coords_2d.y;
      SD.SrcPos(qi, 3) = 0;
      
      SD3D.SrcPos(qi, 1) = layout.docks(didx).optodes(oidx).coords_3d.x;
      SD3D.SrcPos(qi, 2) = layout.docks(didx).optodes(oidx).coords_3d.y;
      SD3D.SrcPos(qi, 3) = layout.docks(didx).optodes(oidx).coords_3d.z;
  end
  
end

% Add detector entries
for mi = 1:SD.nDets
  
  nidx = gldet(mi).node_idx;
  oidx = gldet(mi).optode_idx;
  nid = enum.groups(gi).nodes(nidx).id;
  didx = layout.dockmap(nid);
  
  switch sdstyle
    case 'standard'
      SD.DetPos(mi, 1) = layout.docks(didx).optodes(oidx).coords_3d.x;
      SD.DetPos(mi, 2) = layout.docks(didx).optodes(oidx).coords_3d.y;
      SD.DetPos(mi, 3) = layout.docks(didx).optodes(oidx).coords_3d.z;
      
    case 'flat'
      SD.DetPos(mi, 1) = layout.docks(didx).optodes(oidx).coords_2d.x;
      SD.DetPos(mi, 2) = layout.docks(didx).optodes(oidx).coords_2d.y;
      SD.DetPos(mi, 3) = 0;
      
      SD3D.DetPos(mi, 1) = layout.docks(didx).optodes(oidx).coords_3d.x;
      SD3D.DetPos(mi, 2) = layout.docks(didx).optodes(oidx).coords_3d.y;
      SD3D.DetPos(mi, 3) = layout.docks(didx).optodes(oidx).coords_3d.z;
  end
  
end

% Build MeasList
%
SD.MeasList = zeros(size(glch,1), 4);
SD.MeasList(:,1) = glch(:,1);           % Global source index
SD.MeasList(:,2) = glch(:,3);           % Global detector index
SD.MeasList(:,3) = 1;                   % Always one
SD.MeasList(:,4) = glch(:,2);           % Global wavelength index

% Build event (stimulus) maitrx
%
% LUMO records characters and strings as event markers, which differs from the way that
% stimuli are encoded in NIRS. There is no unambiguous mapping so we implement the approach
% taken in the DOT-HUB toolbox. The following implementation is derived from this code,
% under the GPLv3.
%
% License: https://github.com/DOT-HUB/DOT-HUB_toolbox/blob/master/LICENSE
%
if ~isempty(events)
  
  for i = 1:length(events)
    timeStamp(i) = events(i).timestamp;
    eventStr{i} = events(i).mark;
  end
  
  % Find unique events, maintain order in which they occured
  [tmp, ~, occuranceInd] = unique(eventStr, 'stable');
  
  % Extract condition names
  CondNames = arrayfun(@cellstr, tmp);
  nc = size(CondNames,2);
  
  s = zeros(size(t,1), nc);
  
  for i = 1:nc
    timeStampTmp = timeStamp(occuranceInd == i); %Can have multiple entries
    [~, indTmp] = min(abs(repmat(t, 1, length(timeStampTmp)) ...
      - repmat(timeStampTmp, size(t, 1), 1)));
    s(indTmp,i) = 1;
  end
  
else
  
  
  % No events recorded
  s = [];
  CondNames = {};
  
end

% Sort measurement list and data by (wavelength, source, detector) and apply to data,
% transposing to specified [nt x nc]
[SD.MeasList, perm] = sortrows(SD.MeasList,[4,1,2]);
d = data(gi).chn_dat(perm,:).';

% Add MeasListAct with all channels enabled
SD.MeasListAct = ones(size(SD.MeasList, 1), 1);

% Add MeasListSat with saturation flags
%
% This is a DOT-HUB extension which encodes if a channel is -ever- saturated. The full data
% is exposed int the lumoext structure.
SD.MeasListActSat = data(gi).chn_sat;
SD.MeasListActSat = SD.MeasListActSat(perm,:);

% Apply to 3D structure
if strcmp(sdstyle,'flat')
  SD3D.MeasList = SD.MeasList;
  SD3D.MeasListAct = SD.MeasListAct;
  SD3D.MeasListActSat = SD.MeasListActSat;
  SD3D.SrcPowers = SD.SrcPowers;
end

% A duplicate entry is specified for unknown reasons
ml = SD.MeasList;

% Build output structure
nirs.t = t;
nirs.d = d;
nirs.SD = SD;
nirs.ml = ml;
nirs.s = s;
nirs.CondNames = CondNames;

if strcmp(sdstyle, 'flat')
  nirs.SD3D = SD3D;
end

% Add lumo extension variables
nirs.lumo.chn_sort_perm = perm;

end

