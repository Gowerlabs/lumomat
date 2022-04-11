classdef LumoData
  %LUMODATA Encapsulate the data of a LUMO recording sessions
  %
  % The LUMODATA class implements a high-level API for exploring the enumeration and data of
  % a LUMO recording.
  %
  % Notes:
  %
  %   - A LumoData object can only be created for datasets which can be spectroscopically
  %     indexed, that is to say that all source positions form channels over the same set of
  %     wavelengths.
  %
  %   - If a file is loaded containing data from multiple groups (e.g. a hyperscanning
  %     experiment), a group index must be provided to each function for which data is to be
  %     returned.
  %
  %
  %   (C) Gowerlabs Ltd., 2022
  %
  
  properties (SetAccess = private)
    
    % Canonical data structures
    %
    enum  % Device enumeration
    data  % Raw recording data (incl. ancillary)
    evts  % Events
    
  end
  
  
  properties (Access = private)
    
    % Object-level derived quantities
    %
    ngroups       % Number of groups in the data structure
    
    % Group-level derviced quanties (ngroup x 1)
    %
    nsrcs         % Number of sources
    ndets         % Number of detectors
    nsrcpos       % Number of optodes containing sources
    ndetpos       % Number of optodes containing detectors
    nwls          % Number of wavelengths in the system
    
  end
  
  
  methods
    
    % Constructor
    %
    function obj = LumoData(fn, varargin)
      %LUMODATA Read lumo data from disk
      %
      % [ldat] = LUMODATA(filename, ...)
      %
      % LUMODATA reads a LUMO file from disk, returning an object which encapsulates the
      % system enumeration, layout, channel data, and event markers.
      %
      %   Paramters:
      %
      %   'filename':       The path of the LUMO file to load.
      %
      %   Optional Parameters:
      %
      %   'layout':         When a layout is specified, the embedded layout in the specified
      %                     LUMO file, if present, is ignored, and the alternative
      %                     specification applied. An alternative layout can be specified as
      %                     either:
      %
      %                     string: the filename of a valid lumo layout file, in JSON format
      %                     struct: a layout structure in the format returned by
      %                             lumofile.read_layout (see function help for details)
      %
      %                     If the provided layout has been constructed by the user, for
      %                     example, based upon measurements of a physical layout, entries
      %                     must be present for each occupied dock in the recording with the
      %                     apprporiate dock ID.
      %
      %   'ignore_memory':  A logical value which indicates if the function should ignore 
      %                     any potential performance issues arising from loading large
      %                     amounts of data. Defaults to false.
      %
      %
      %   (C) Gowerlabs Ltd., 2022
      %
      
      [~,~,ext] = fileparts(fn);
      
      switch lower(ext)
        case '.lumo'
          [obj.enum, obj.data, obj.evts] = lumofile.read_lumo(fn, varargin{:});
        case 'lufr'
          error('LUFR file loading is not implemented');
          % [obj.enum, obj.data, obj.evts] = lumofile.read_lufr(fn, varargin{:});
        otherwise
          error('Unknown file extension %s', ext');
      end
      
      warning('REQUIRE LAYOUT, REQUIRE SPECTROSCOPIC, COMPLETE DERIVED');
      
    end
    
    % File output
    %
    function nirs = write_NIRS(obj, fn, varargin)
      % WRITE_NIRS Covert LUMO data to NIRS format and write to disk
      %
      %   [nirs] = WRITE_NIRS(fn)
      %
      % See LUMOFILE.write_NIRS for details.
      %
      [nirs] = lumofile.write_NIRS(fn, obj.enum, obj.data, obj.evts, varargin{:});
    end
    
    function write_SNIRF(obj, fn, varargin)
      % WRITE_NIRS Covert LUMO data to NIRS format and write to disk
      %
      %   [nirs] = WRITE_NIRS(filename)
      %
      % See LUMOFILE.write_SNIRF for details.
      %
      lumofile.write_SNIRF(fn, obj.enum, obj.data, obj.evts, varargin{:});
    end
    
    % Channel query
    %
    function [info] = chn_info(obj, cidx, varargin)
      % CHN_INFO Return information on a particular channel, by index
      %
      %   [info] = chn_info(index, ...)
      %
      % CHN_INFO returns a complete description of a channel by its index. The information
      % structure contains the following fields:
      %
      %   info.idx:             Channel index
      %       .src_node_id:     ID of the node containing the source of this channel
      %       .src_idx:         Index of the source within the node
      %       .src_optode_name: Name of the source optode(e.g. 'A', 'B', 'C')
      %       .src_wl:          Wavelength of the source [nm]
      %       .src_coord_2d:   (*) Two-dimensional co-ordinates of the source
      %       .src_coord_3d:   (*) Three-dimensional co-ordinates of the source
      %       .det_node_id:     ID of the node containing the detector of this channel
      %       .det_idx:         Index of the detector within the node
      %       .det_optode_name: Name of the detector optode (e.g. '1', '2', '3', '4')
      %       .det_coord_2d:   (*) Two-dimensional co-ordinates of the source
      %       .det_coord_3d:   (*) Three-dimensional co-ordinates of the source
      %
      % (*) Optional fields which may be empty if the requisite information is unavailable.
      %
      % Notes:
      %
      %  - Multiple indices can be provided, resultig in a an output which is an array of
      %  structures.
      %
      %  - If the data object contains multiple groups, the argument pair ('group', index)
      %  will be required to determine the correct enumeration.
      %
      %
      %   (C) Gowerlabs Ltd., 2022
      %
      
      % Group
      gidx = method_group_idx(obj, varargin{:});
      
      % Over each channel
      for i = 1:length(cidx)
        
        % Channel
        ch = obj.enum.groups(gidx).channels(cidx(i));
        
        % Source
        src_node_idx = ch.src_node_idx;
        src_node = obj.enum.groups(gidx).nodes(src_node_idx);
        src_idx = ch.src_idx;
        src_optode_idx = src_node.srcs(src_idx).optode_idx;
        src_wl = src_node.srcs(src_idx).wl;
        src_dock_idx = obj.enum.groups.layout.dockmap(src_node.id);
        src_dock = obj.enum.groups.layout.docks(src_dock_idx);
        src_optode = src_dock.optodes(src_optode_idx);

        % Detector
        det_node_idx = ch.det_node_idx;
        det_node = obj.enum.groups(gidx).nodes(det_node_idx);
        det_idx = ch.det_idx;
        det_optode_idx = det_node.dets(det_idx).optode_idx;
        det_dock_idx = obj.enum.groups.layout.dockmap(det_node.id);
        det_dock = obj.enum.groups.layout.docks(det_dock_idx);
        det_optode = det_dock.optodes(det_optode_idx);

        % Build the information
        info(i) = struct('idx', cidx(i),...
          'src_node_id', src_node.id, ...
          'src_idx', src_idx, ...
          'src_optode_name', src_optode.name, ...
          'src_wl', src_wl, ...
          'src_coord_2d', [src_optode.coord_2d.x src_optode.coord_2d.y], ...
          'src_coord_3d', [src_optode.coord_3d.x src_optode.coord_3d.y src_optode.coord_3d.z], ...
          'det_node_id', det_node.id, ...
          'det_idx', det_idx, ...
          'det_optode_name', det_optode.name, ...
          'det_coord_2d', [det_optode.coord_2d.x det_optode.coord_2d.y], ...
          'det_coord_3d', [det_optode.coord_3d.x det_optode.coord_3d.y det_optode.coord_3d.z]);
        
      end

    end
    
    function [idx] = chn_find(obj, varargin)
      % CHN_FIND Find channel indices by properties
      %
      %   [indices] = chn_find('property', value, ...)
      %
      % CHN_FIND returns a list of indices of channels which match all of the specified
      % properties. Each property is is (name, value) pair.
      % 
      % Supported filters:
      %
      %   'src_node_id':      ID of the source node
      %   'det_node_id':      ID of the detector node
      %   'intratile':        When true, filters channels to those within a tile
      %   'intertile':        When true, filters channels to those across tiles
      %   'wavelength':       Filters by wavelength (specified in nm)
      %
      % Example:
      %
      % Find all channels of wavelength 735nm, where the source is on node 1, and the
      % detector is on node 8:
      %
      % >> ld.chn_find('wavelength', 735, 'src_node_id', 1, 'det_node_id', 8).'
      % 
      % ans =
      % 
      %     29
      %     30
      %     31
      %     32
      %    125
      %    126
      %    127
      %    128
      %    221
      %    222
      %    223
      %    224

      %
      %
      % Notes:
      %
      %  - If the data object contains multiple groups, the argument pair ('group', index)
      %  will be required to determine the correct enumeration.
      %
      %
      %   (C) Gowerlabs Ltd., 2022
      %
      
      %%% TODO
      %
      %   'src_name':         Name of the source (e.g., 'A', 'B', 'C')
      %   'det_name':         Name of the detector (e.g., '1', '2', '3', '4')
         
      % Group
      gidx = method_group_idx(obj, varargin{:});
      
      % Parse parameters
      p = inputParser;
      addOptional(p, 'src_node_id', [], @isnumeric);
      addOptional(p, 'det_node_id', [], @isnumeric);
      addOptional(p, 'src_name', [], @ischar);
      addOptional(p, 'det_name', [], @ischar);
      addOptional(p, 'intratile', false, @islogical);
      addOptional(p, 'intertile', false, @islogical);
      addOptional(p, 'wavelength', [], @isnumeric);
      parse(p, varargin{:});
      
      f_src_node_id = p.Results.src_node_id;
      f_det_node_id = p.Results.det_node_id;
%       f_src_name = p.Results.src_name;
%       f_det_name = p.Results.det_name;
      f_intratile = p.Results.intratile;
      f_intertile = p.Results.intertile;
      f_wavelength = p.Results.wavelength;
      
      % Prepare variables
      channels = obj.enum.groups(gidx).channels;
      nodes = obj.enum.groups(gidx).nodes;
      
      % Start with the full set
      cidx = true(length(channels), 1).';
            
      % Apply filters...
      if ~isempty(f_src_node_id)
        src_node_ids = [nodes([channels.src_node_idx]).id];
        if isscalar(f_src_node_id)
          fidx = (src_node_ids == f_src_node_id);
        else
          fidx = any(src_node_ids == f_src_node_id);
        end
        cidx = cidx & fidx;
      end

      if ~isempty(f_det_node_id)
        det_node_ids = [nodes([channels.det_node_idx]).id];
        if isscalar(f_det_node_id)
          fidx = (det_node_ids == f_det_node_id);
        else
          fidx = any(det_node_ids == f_det_node_id);
        end
        cidx = cidx & fidx;
      end
      
      if f_intratile
        fidx = [channels.src_node_idx] == [channels.det_node_idx];
        cidx = cidx & fidx;
      end
            
      if f_intertile
        fidx = [channels.src_node_idx] ~= [channels.det_node_idx];
        cidx = cidx & fidx;
      end
      
      if ~isempty(f_wavelength)
        
        wls = zeros(length(channels),1).';
        for i = 1:length(channels)
          wls(i) = nodes(channels(i).src_node_idx).srcs(channels(i).src_idx).wl;
        end
        if isscalar(f_wavelength)
          fidx = wls == f_wavelength;
        else
          fidx = any(wls == f_wavelength.');
        end
        cidx = cidx & fidx;        
      end
      
      idx = find(cidx);
            
    end
    
    function [data] = chn_data(obj, idx, varargin)
      % CHN_DATA Return intensity time series for specified channel indices
      %
      %   [data] = chn_data(indices, ...)
      %
      % CHN_DATA returns an [nc x nt] array of channel data over the [nc] channel indices
      % over all [nt] time points.
      %
      %  - If the data object contains multiple groups, the argument pair ('group', index)
      %  will be required to determine the correct enumeration.
      %
      %
      %   (C) Gowerlabs Ltd., 2022
      %
      
      % Group
      gidx = method_group_idx(obj, varargin{:});
      
      % Return data
      data = obj.data(gidx).chn_dat(idx, :);
            
    end
    
    function [t] = chn_time(obj)
      % CHN_TIME Return time vector for channel data 
      %
      %   [t] = chn_time()
      %
      % CHN_TIME returns an [1 x nt] vector of time points at which each intensity measure
      % was recorded in units of seconds.
      %
      %   (C) Gowerlabs Ltd., 2022
      %
      
      % Group (all groups will be at the same rate)
      gidx = 1;
      
      % Return data
      t = (0:(obj.data(gidx).nframes - 1))*obj.data(gidx).chn_dt*1e-3;
      
    end     
      
    
    % Reindexing
    %
    function [ch, sp, dp, wl] = flatten(obj, varargin)
      
      % Group
      gidx = method_group_idx(obj, varargin{:});
      
      % Global spectroscopic
      % ch -> (source, detector, wavelength)
      
      % Global
      % ch -> (source, detector)
      
      % Local spectroscopic
      % ch -> (source, detector)
      
      ch=[];
      sp=[];
      dp=[];
      wl=[];
    end
    
  end
  
  methods (Access = private)
    
    % Helper to manage group indexing
    function [gidx] = method_group_idx(obj, varargin)
      
      if obj.ngroups > 1
        
        % Require specification of group index
        if ~strcmpi(varargin{1}, 'group')
          error(['LumoData object contains information from more than one group. Specify' ...
            'the required group using the argument pair (''group'', group_idx)']);
        else
          gidx = varargin{2};
          if (gidx < 1) || (gidx > length(obj.enum.groups))
            error('Group index %d out of range', gidx);
          end
        end
        
      else
        
        % Default
        gidx = 1;
        
      end
      
    end
    
  end
  
  
  
end


