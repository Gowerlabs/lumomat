# lumomat

MATLAB tools for LUMO data.

`lumomat` allows one to read LUMO output files, explore their contents, and export to common file formats. 

# Installation

1. Download the latest version from the [release page](https://github.com/Gowerlabs/lumomat/releases)
2. Unzip in your favourite directory
3. Add the  `lumomat` directory to your MATLAB path

# Quickstart

Load a `.lumo` or a a `.lufr` file:

```
data = LumoData(filename)
```

Convert the file to a SNIRF file:

```
data.write_SNIRF(filename);
```

Get information about channel 342 of your data:

```
>> info = data.chn_info(342)

info = 

  struct with fields:

                idx: 342
        src_node_id: 2
            src_idx: 4
    src_optode_name: 'A'
             src_wl: 850
      src_coords_2d: [94.3617 75.6020]
      src_coords_3d: [87.3064 179.5761 110.9974]
        det_node_id: 2
            det_idx: 2
    det_optode_name: '2'
      det_coords_2d: [102.0520 81.9720]
      det_coords_3d: [79.4308 175.0386 115.1332]
```

Get information about the channel from node ID 4, source optode 'B', wavelength 850nm, to node ID 7, detector optode '2' (note that the optode naming follows the tile description in the LUMO user manual):

```
>> data.chn_info(data.chn_find('src_node_id', 4, 'src_optode_name', 'B', 'wavelength', 850, 'det_node_id', 7, 'det_optode_name', '2'))


ans = 

  struct with fields:

                idx: 1035
        src_node_id: 4
            src_idx: 5
    src_optode_name: 'B'
             src_wl: 850
          src_power: 80
      src_coords_2d: [41.8617 45.3010]
      src_coords_3d: [135.4298 179.5515 76.8679]
        det_node_id: 7
            det_idx: 3
    det_optode_name: '3'
      det_coords_2d: [146.8617 56.1110]
      det_coords_3d: [41.0043 173.0580 89.2694]

```

Find all of the channel indices where the source node ID is 1, the detector node ID is 8, and the wavelength is 850nm:

```
>> idx = data.chn_find('wavelength', 850, 'src_node_id', 1, 'det_node_id', 8).'

ans =

    77
    78
    79
    80
   173
   174
   175
   176

   269
   270
   271
   272
```

Find the name of the source optode with the first result from our previous query:

```
>> data.chn_info(idx(1)).src_optode_name

ans =

    'A'
```

And the channel data for this channel:

```
>> data.chn_data(idx(1)).'

ans =

  978×1 single column vector

   1.0e-04 *

    0.3666
    0.6109
    0.2772
    0.4858
    .
    .
    .
```

# Introduction

LUMO is a high-density, wearable, and modular instrument for diffuse optical tomography (DOT) and (functional) near-infrared spectroscopy ((f)NIRS). 

In principle, the data produced in a DOT experiment is straightforward: a single measurement of intensity is made across every channel of the system many times per second, resulting in a big matrix of data. In practice, there is a significant amount of book-keeping required in order to make sense of the relationship between that big matrix of numbers and the physical geometry of the system.

This is particularly true for LUMO, owing to its modular nature. The purpose of this package is to:
 - normalise various LUMO file formats, and file format versions to a common representation suitable for interpretation and analysis,
 - allow export of the data in file formats used by common analysis packages.
  

# Nomenclature

Various terms are used throughout this guide and in the package itself. Many will be familiar to practitioners in (f)NIRS and DOT, others are specific to LUMO, and the most dangerous are ambiguous.

 - *Tile* / *Node* / *Module*: a LUMO tile is a single optoelectronic module containing a set of optical sources, detectors, motion processing units, and so on. Lumo tiles are plugged into docks to form a DOT system. In the logical representation of the system, a tile is called a 'node' (all tiles are nods, but not all nodes are tiles, but for most practical situations they are synonymous). A node is the same concept as a 'module' in the parlance of the SNIRF file format.
 - *Dock*: a dock is the receptacle into which a LUMO tile is installed. Multiple docks are normally connected together in a piece of headgear, or a patch, which is applied to the subject of the experiment.
 - *Group* / *Cap*: most fNIRS and DOT experiments involve measurements of the haemodynamics of the brain, and thus most systems are formed of headgear called caps into which multiple docks are installed. But there are many other interesting targets which can be measured, so the generic name for a group of docks which form a single measurement system is, imaginatively, a 'group'.
 - *Optode*: an optode is a single location from which light is emitted by one or more sources, and/or detected by one or more detectors. Conceptually there is a one-to-one mapping between the optode on a tile, and the optode on a dock: this is the link between logical sources and physical locations.
 - *Source*: LUMO considers a source to mean a single optical emitter at a specific wavelength. There are usually multiple sources at different wavelengths co-located at a single optode in order to enable spectroscopy. *This meaning may differ from some file formats and analysis packages which think of a source as an optode, with multiple wavelengths.*
 - *Detector*: a detector is a device which measures optical intensity. Like sources, detectors are associated with a particular optode.
 - *Channel*: a channel is formed from a single source-detector pair.
 - *Layout*: a layout describes the docks present in a group, and each of the optodes which belong to the dock. When a cap is built, a template layout is created which provides the positions of all docks in a group as measured on a suitable phantom. Subject specific layouts might be measured by the user (and may consist of a subset of the docks in the group).
 - *Frame*: a frame is a measurement of all active channels made at a single point in time.

# Using lumomat

The main interface to `lumomat` is the `LumoData` class, which is constructed from a LUMO output file:

```
data = LumoData(filename)
```

A `LumoData` object has a number of methods which permit exploration of the data contained in the file. We summarise the methods here, see online help (or the Quickstart section) for further details.

Channels can be explored:

 - `chn_find`: locate channels by filtering on parameters
 - `chn_info`: collate and return all information regarding a particular channel
 - `chn_data`: get the raw data from a particular channel

The system enumeration can be re-indexed:

 - `reindex_global`: convert system description to global spectroscopic format, used by many DOT and (f)NIRS analysis programs, see the Indexing section below for further details.

The data can be exported (with LUMO extensions, see below):

 - `write_SNIRF`: output to the .SNIRF format 
 - `write_NIRS`: output to the .NIRS format

## Layouts

The physical layout of a LUMO system is determined from a layout file which contains the position of each optode in each dock of a given group, and optional physiological landmarks. During manufacture, these positions are measured on a suitable phantom and a *default* layout file is produced. This layout file is copied into the output file when recording.

Note: some versions of the LUMOview software do not embed a layout file in the `.lumo` output file as described above. If no layout file is embedded, a warning will be printed. In this case:
-  the `lumofile.merge_layout` function can be used to embed a layout file in a copy of a `.lumo` recoding. When this function is called without specifying a layout file, it will search the standard installation directories for the layout file included with the LUMOview software
 -  the layout file may be specified as described for the provision of subject-specific layouts
  
It is often the case that for the purposes of data analysis, practitioners will wish to provide subject-specific layout information measured during an experiment. This information can be used instead of the default layout:

```
>> data = LumoData(filename, 'layout', custom_layout);
```

where `custom_layout` can be one of:
 - a structure matching that described in the layout section of the low-level API description, or, 
 - an JSON file following the syntax of the default layout file (contact Gowerlabs for details)

For a custom layout to be usable it must contain the locations of all docks for which nodes were present. This may be a subset of the complete set of docks in a group, in the case that some are not occupied.

All data export methods write additional metadata to, e.g., NIRS or SNIRF files which permit identification of the node and optode corresponding to each source or detector position. As such, it is also possible for the user to modify the positions using measured data after export.

## Indexing

The canonical enumeration of a LUMO system uses (node-) local indexing, and the sources are not inherently combined by wavelength to form optodes. As such, channels are represented internally by a 4-tuple:

*ch = (source node, source index, detector node, detector index)*

This representation permits complete flexibility in the description of channels in the system, and maintains the link between channels, nodes, and docks, which is useful in a modular system. 

Many data formats and analysis tools for (f)NIRS and DOT use an alternative format to index the channels of the system. The most common choice is the global spectroscopic scheme, in which the channels are indexed by their optodes and wavelength. In this scheme, a channel in defined by a 3-tuple:

*ch = (source optode, source wavelength, detector optode)*

This indexing scheme is less flexible, but it is sane because:

 - most experiments involve spectroscopy, so it is a fair assumption that all source optodes will have associated with them sources at every wavelength
 - most systems do not collocate sources and detectors in the same optode, so it is reasonable to assign an optode to be either belonging to sources, or detectors

The canonical enumeration of a LUMO system can usually be converted to global spectroscopic indexing, and `lumomat` enforces that this is the case when you use the `LumoData` API. A method is provided to re-index the enumeration, returning the information required for further analysis or export to custom file formats:

```
[global_chns, src_optodes, det_optodes, global_wls] = data.reindex_global();
```

The `global_chns` output is an array of channel descriptors:

```
>> global_chns

global_chns = 

  1×3456 struct array with fields:

    src_optode_idx
    det_optode_idx
    wl_idx
```

Where the `src_optode_idx`, `det_optode_idx`, and `wl_idx` fields index respectively into the remaining outputs. For example, we can inspect channel 932:

```
>> src_optodes(global_chns(932).src_optode_idx)

ans = 

  struct with fields:

      node_id: 4
         name: 'A'
    coords_2d: [32.5000 61.5160]
    coords_3d: [138.9866 164.5234 87.4543]

>> global_wls(global_chns(932).wl_idx)

ans =

   850

>> det_optodes(global_chns(932).det_optode_idx)

ans = 

  struct with fields:

      node_id: 5
         name: '3'
    coords_2d: [76.8617 64.9910]
    coords_3d: [104.2389 185.0212 102.2458]
```

The ordering of the channels is maintained, so we can compare with the internal canonical enumeration:

```
>> data.chn_info(932)

ans = 

  struct with fields:

                idx: 932
        src_node_id: 4
            src_idx: 4
    src_optode_name: 'A'
             src_wl: 850
       src_coords_2d: [32.5000 61.5160]
       src_coords_3d: [138.9866 164.5234 87.4543]
        det_node_id: 5
            det_idx: 4
    det_optode_name: '3'
       det_coords_2d: [76.8617 64.9910]
       det_coords_3d: [104.2389 185.0212 102.2458]

```       

## Channel saturation and frame errors

In common with most systems, it is possible for channels to become saturated owing to subject movement following the setting of optical source powers. The multiplexing strategy employed by LUMO is such that saturated channels may appear to have a *lower* signal than when unsaturated. To determine if a channel is saturated, consult the appropriate entry in the `chn_sat` field of the `LumoData` object. 

The format of the channel saturation data in the `LumoData` output depends on the version of the software which was used to record the data. Recent versions of the LUMO software will provide a saturation indication for every channel, for every frame, in which case the `chn_sat` field will have the dimensions `<no. frame x no. channels>`. Older versions of the software may not output the channel saturation data at all or only indicate if a channel is saturated at any time during the recording.

In rare circumstances, high CPU or USB load on the acquisition computer may interrupt communication with LUMO. Recent versions of the LUMO software will tolerate small interruptions without terminating a recording, but the resultant frames should be excluded from analysis. To determine if a frame should be excluded, check for a non-zero entry in the `err_cnt` field of the `LumoData` object.

## Auxiliary data

The `LumoData` object may contain additional auxiliary data, such as measurements from integrated accelerometers, gyroscopes and temperature sensors:

```
>> data

ans = 

  struct with fields:

         chn_dat: [24×1503 single]
          chn_dt: 80
         chn_fps: 12.5000
         chn_sat: [24×1503 logical]
         err_cnt: [1503×1 double]
         nframes: 1503
           nchns: 24
       node_temp: [1×1503 single]
     node_mpu_dt: 0.0100
    node_mpu_fps: 100
        node_acc: [1×3×12024 double]
        node_gyr: [1×3×12024 double]
```

When present, the accelerometer and gyroscope data has dimensions of `<no. tiles x 3 x no. time>` where the second dimension is indexed over the x, y, and z-axes. The time vector for such motion data can be computed using the `node_mpu_dt` field.

## SNIRF output

The [SNIRF file format](https://github.com/fNIRS/snirf) is a recent specification supported directly by a number of modern analysis tools such as [MNE-NIRS](https://github.com/mne-tools/mne-nirs), [Fieldtrip](https://www.fieldtriptoolbox.org/), [Homer3](https://github.com/BUNPC/Homer3), and [NIRS-toolbox](https://github.com/huppertt/nirs-toolbox).

SNIRF uses HDF5 as its underlying data storage method. HDF5 is a mature and stable file format which has been proven in numerous large scale experimental systems. The format has excellent support across multiple languages and on various platforms. *For these reason, we recommend the use of SNIRF as an archival format for LUMO data*.

Whilst SNIRF has some support for local indexing (similar to the canonical indexing used internally), this may not be the case for tools which use the format. As such, we export to SNIRF using global spectroscopic indexing, as described previously. We also include metadata to permit easy mapping of the channels of the system back to a local format, such that it is straightforward to determine, e.g., the nodes involved in a given channel.

To construct a SNIRF structure and write a SNIRF file:

```
snirf = data.write_SNIRF(filename);
```

The returned structure represents the contents of the SNIRF output file as a nested MATLAB structure, where indexed groups are replaced with arrays of structures. This data structure can be used directly if a globally indexed spectroscopic representation is of use.

We can examine the root of the first measurement group `nirs{1}`:

```
>> snirf.nirs(1)

ans = 

  struct with fields:

    metaDataTags: [1×1 struct]
            data: [1×1 struct]
           probe: [1×1 struct]
            stim: [1×4 struct]
```

Look at the metadata:

```
>> snirf.nirs(1).metaDataTags

ans = 

  struct with fields:

           SubjectID: 'Subject Unknown'
     MeasurementDate: 'unknown'
     MeasurementTime: 'unknown'
          LengthUnit: 'mm'
            TimeUnit: 'ms'
       FrequencyUnit: 'Hz'
     sourcePowerUnit: 'percent'
    ManufacturerName: 'Gowerlabs'
                Mode: 'LUMO'
                lumo: [1×1 struct]
```

Or examine the nature of the measurement list structure:

```
>> snirf.nirs(1).data(1).measurementList

ans = 

  1×384 struct array with fields:

    sourceIndex
    detectorIndex
    wavelengthIndex
    dataType
    dataTypeIndex
    sourcePower
```

### LUMO conventions

Some optional output fields are formatted in a manner specific to the LUMO system:

 - `/nirs{i}/probe/sourceLabels`: source labels are formatted as "N\<node ID\>-\<name\>\<wavelength\>", for example, source A, wavelength 850, on node ID 7 will have the name "N7-A850".
 - `/nirs{i}/probe/detectorLabels`: detector labels are formatted as "N\<node ID\>-\<name\>", for example, detector 3, on node ID 28 will have the name "N28-3".

### LUMO stimulus data

The SNIRF format expects stimuli to be be recorded as set of conditions, where each condition can have multiple trials which are recorded as a start time, a duration, and an amplitude. 

The LUMO system does not explicitly record stimuli in this format, instead it provides for the recording of *event markers*, each of which consists of a string (the name of the event) and the time of occurrence. In a practical (f)NIRS experiment event markers may be used with a-priori knowledge of the stimulus duration, or different markers might be used to indicate the start and end of a stimulus condition.

To encode event markers in the SNIRF stimulus format, all events which share the same marker are recorded as individual stimulus conditions. Every incidence of an event marker is recorded as a stimulus trial with the appropriate start time, a zero duration, and a unit amplitude. *The user is responsible for modifying the resultant data according to the specific experimental paradigm.*

### LUMO additional metadata

This package writes additional fields in the SNIRF metadata, as permitted by the SNIRF specification. The following fields are located under `nirs(i)/metaDataTags/`:

 - `lumomatVersion`: the version of the `lumomat` package which wrote the SNIRF file, a string representation of a semantic version number.
 - `saturationFlags`: a vector of integers in which a non-zero value in the `i`th element indicates that saturation of the `i`th channel occurred at some time during the recording. Transient saturation can occur during, e.g., movement, so this global flag can exclude many channels which are viable for the vast majority of the recording. Some versions of the LUMO software will export a time-series of saturation flags (see auxiliary measurements) to enable more granular channel filtering.
 - `errorFlags`: a vector of integers in which a non-zero value in the `i`th element indicates that a data reception error occurred in the `i`th frame of the recording. This frame will likely contain NaN entries for one or more data types, from one or more tiles.
 - `groupName`: the name of group upon which the data was acquired (e.g. the cap serial number).

Additional extended metadata can also be written to the output file, however the resultant file is not compliant with the SNIRF specification. To enable writing the extended metadata:

```
snirf = data.write_SNIRF(filename, 'meta', 'extended');
```

This will result in the following fields being written under the  `nirs(i)/metaDataTags/lumo/` group:

 - `canonicalMap`: is a flattened and abbreviated version of the canonical enumeration used internally, the values of this matrix can be used to index into the `nodes` and `docks` indexed groups. The values can also be used on their own in order to restore locality information (e.g. the nodes to which channels belong) from the global enumeration stored in the `probe` group. The rows of the matrix are each 1-based indices:
   1. source node index, used to index the `nodes{i}` structure
   2. source dock index, used to index the `docks{i}` structure
   3. source optode index, used to index the `docks{i}
   4. /optodePosXD` arrays
   5. source wavelength, in nm
   6. detector node index
   7. detector dock index 
   8. detector optode index
 - `nodes{i}`: an indexed group of information about nodes, such as their `id`
 - `docks{i}`: an indexed group of all docks in the layout. Note that unlike the global enumeration stored in the `probe` group, which only contains information about the optode positions which are in use, the dock structure contains information about all the docks in the system. Each dock contains a number of datasets, including:
   - information such as the optode `name`
   - an array of optode locations in flattened two-dimensional and three-dimensional co-ordinates `optodePosXD`


### LUMO auxiliary measurements

Auxiliary measurements such as motion data and temperature are also written to the output file, when available. The format of this data is not yet stable as it depends upon the resolution of a [query](https://github.com/fNIRS/snirf/issues/107) regarding the SNIRF specification.

Each auxiliary measurement is located in an individual `/nirs(i)/aux(j)` field. The `/nirs{i}/aux{j}/name` field will be set as decried to identify the data.

  - (*) Temperature (`temperature`): matrix of dimensions `<time x nodes>` containing the internal tile temperature for each node over the course of the experiment. 
  - (*) Saturation (`saturationFlags`): matrix of dimensions `<time x channels>` where a non-zero value indicates that the channel was saturated at the selected time point. This is a more detailed version of the channel saturation summary written into the LUMO metadata.
  - (*) Motion data (`accel_x/y/z`, `gyro_x/y/z`): a set of six matrices each of dimensions `<time x nodes>`, with names according to the SNIRF specification. Acceleration is stored in units of `g`, and gyroscope data in `deg/s`. Note that motion data is not unfiltered, and post-processing will be required.
  
*All auxiliary fields are optional and only recorded by more recent software versions of LUMO. To update your software in order to capture this data, contact Gowerlabs.*


### Notes

 - When a layout field contains physiological landmarks, these will be stored in the appropriate `landmarkLabels` and `landmarkPos3D` datasets, but `landmarkPos2D` is optional and may not be present.


## NIRS output

The NIRS file format is based upon MATLAB output files, and was originally used in the HOMER2 software. Owing to its legacy, a number of existing processing pipelines might most easily be used with data in this format. Whilst a specification is provided for valid NIRS files, variations on the format are commonly encountered. We recommend that when possible, users prefer the SNIRF format for their data archival and analysis.

To construct a NIRS structure and write a NIRS file:

```
nirs = data.write_NIRS(filename);
```

The NIRS file contains an `SD` structure which describes the physical configuration of the sources and detectors of the system. In some references this layout is considered to be two-dimensional, though the specification provides for three-dimensional co-ordinates. The `write_NIRS` method provides an option which allows the `SD` structure to be written in two different formats to accommodate varying use cases:

 - `standard` (default): the SD structure is built using the proper three-dimensional co-ordinates of the system
 - `flat`: the SD structure is built using a flattened two-dimensional layout, and an additional `SD3D` structure is included in the NIRS file which contains the three-dimensional data. Files written in this format are compatible with the DOT-HUB toolbox.

 To select a scheme, pass the appropriate argument pair, e.g.: 
 
 ```
 data.write_NIRS(filename, 'sd_style', 'flat');
```

*Note: The channel list in a NIRS file is re-indexed such that it is sorted by (wavelength, source index, detector index), as this is assumed in some analysis software.*

### Additional LUMO specific fields:

 - `SD.MeasListActSat`: a logical vector (or matrix) indicating if a channel is saturated. If this field is a vector, a non-zero (or logical true) entry in the `i`th index indicates that the `i`th channel was saturated at some point in the recording. Since transient effects such as movement can cause temporary saturation, the use of a single flag for each channel can cause the loss of a number of channels which are useful for a large proportion of the recording. More recent versions of the LUMOview software will output data which permits allows this field to be output as a matrix of values such that saturation can be identified per-channel, per-frame, allowing for more granular exclusion of saturated data. Contact Gowerlabs to update your software if this feature is desired.
 - `SD.SrcPowers`: a matrix of source powers expressed in percent, indexed by the global source index and wavelength. 
 - `ErrorFlags`: a vector of integers in which a non-zero value in the `i`th element indicates that a data reception error occurred in the `i`th frame of the recording. This frame will likely contain NaN entries for one or more data types, from one or more tiles.

# Low-level functional API

In most circumstances, users should choose the high-level object based API in order to access their data. The high-level wrapper is implemented on-top of a low-level functional API, namespaced in MATLAB packages. 

No guarantees are made regarding the stability of the low-level API, so users must take account of versioning when calling the low-level API.

To begin, load a `.lumo` file:

```
>> [enum, data, events] = lumofile.read_lumo('sample.lumo');
Loading LUMO file sample.lumo
LUMO file (sample.lumo): version 0.3.0
Constructing canonical enumeration...
LUMO file (sample.lumo) enumeration contains 12 tiles, 3456 channels
LUMO file (sample.lumo) events file contains 1 entries
LUMO file (sample.lumo) assigned layout contains 12 docks, 84 optodes
LUMO file loaded in 11.0s
```

Whilst loading the file, information will be printed regarding the contents of the file. Some `.lumo` files might not contain an embedded layout file. If this is the case a warning will be printed.

```
Warning: The specified LUMO file (group C002N / 20155487) does not contain an embedded layout file, and no layout has been specified when calling this function. The returned enumeration will lack layout information, and it will not be possible to convert this file to formats which require a layout. Specify an appropriate layout file to suppress this warning.
```

Follow the guidance, or read more about Layouts in a later section.

Similar functions exist to read `.lufr` files:

```
>> [enum, data, events] = lumofile.read_lufr('sample.lufr');
```

## Enumeration, Enumeration, Enumeration

When you load a `.lumo` or a `.lufr` file, `lumomat` constructs a canonical enumeration which provides a full description of the system. Navigating the enumeration allows the data to be interpreted and manipulated into the form required for further analysis.

The file we have loaded was recorded on a system using 12 tiles, resulting in 3456 channels. The enumeration contains all the information required to interpret the raw data. We can examine its contents:

```
>> enum

enum = 

  struct with fields:

       hub: [1×1 struct]
    groups: [1×1 struct]
```

The `hub` field contains information about the system that is not typically used during data analysis. The `groups` field is a structure describing the group, or cap. When hyper-scanning on a single hub (e.g. multiple groups are recorded at once), you must choose the group of interest when loading the data. Hence, `groups` is scalar and no indexing is required:

```
>> enum.groups

ans = 

  struct with fields:

         uid: '0x0000000001338cbe'
       nodes: [1×12 struct]
    channels: [1×3456 struct]
      layout: [1×1 struct]
```

The `uid` field is a unique identifier which determines the group upon which the recording was made, this can be used to identify an appropriate layout, though this is not enforced.

## Nodes

The `nodes` array contains information about each node that formed the system, including their sources and detectors. For example, the eighth entry:

```
>> node_idx = 8;
>> enum.groups.nodes(node_idx)

ans = 

  struct with fields:

          id: 8
    revision: 2
       fwver: '1.1.0'
        srcs: [1×6 struct]
        dets: [1×4 struct]
     optodes: [1×7 struct]
```

This node happens to have an `id` of 8, as well as being located in the eighth element of the array. This is because the system was fully populated. In general, the node ID and the index in the array do not correspond. We can see from the dimensionality of the `srcs` field that this node has six sources. Let's examine source one:

```
>> src_idx = 1;
>> enum.groups.nodes(node_idx).srcs(src_idx)

ans = 

  struct with fields:

            wl: 735
    optode_idx: 5
         power: 80
```

And we can extract information about all of them, for example, the unique wavelengths:

```
>> unique([enum.groups.nodes(node_idx).srcs.wl])

ans =

   735   850
```

We note that this source has an `optode_idx = 5`. Whenever a field is named ending `_idx` this indicates that the value can be used as (1-based) index into another field. For example, we can determine the optode to which this source is associated:

```
>> optode_idx = enum.groups.nodes(node_idx).srcs(src_idx).optode_idx;
>> enum.groups.nodes(node_idx).optodes(optode_idx)

ans = 

  struct with fields:

    name: 'A'
    type: 'S
```

The name of the optode matches that used in the user manual for the LUMO system, it is only used for display purposes. This optode is of type 'S' which means that it is only associated with optical sources. This is of significance when we discuss global indexing.

We can physically locate this source by looking into the layout. Each dock in the layout shares an ID with the node it accommodates. A mapping is provided to convert an ID to an index into the array of docks:

```
>> node_id = enum.groups.nodes(node_idx).id;
>> dock_idx = enum.groups.layout.dockmap(node_id);
```

Once we have its index we can examine the dock, and the particular optode:

```
>> dock = enum.groups.layout.docks(dock_idx);
>> optode = enum.groups.layout.docks(dock_idx).optodes(optode_idx)

optode = 

  struct with fields:

        name: '1'
   coords_2d: [1×1 struct]
   coords_3d: [1×1 struct]
```

And finally determine its location:

```
>> optode.coords_3d

ans = 

  struct with fields:

    x: 32.2093
    y: 177.5146
    z: 56.1895
```

We will talk more about layouts in a later section.

## Channels

Each channel in the system is defined by a single source-detector pair. A particular source or detector is indexed locally by the node to which it belongs. The channels array thus has the following form:

```
>> enum.groups.channels

ans = 

  1×3456 struct array with fields:

    src_node_idx
    src_idx
    det_node_idx
    det_idx
```

Since we have already developed significant affection for the first source on the eighth node, let's pick a channel for which it provides the illumination:

```
>> channels = enum.groups.channels;
>> find(([channels.src_node_idx] == node_idx) & ([channels.src_idx] == src_idx)).'

ans =

        2017
        2018
        2019
        2020
        .
        .
        .
        2061
        2062
        2063
        2064
```

There are quite a lot, let's pick channel 2061:

```
>> channels(2061)

ans = 

  struct with fields:

    src_node_idx: 8
         src_idx: 1
    det_node_idx: 12
         det_idx: 1
```

The source node and source indices should not come as a surprise. The wavelength of this channel is:

```
>> enum.groups.nodes(node_idx).srcs(src_idx).wl

ans =

   735
```

Following the same procedure as before, we can determine the location of the detector:

```
>> det_node = enum.groups.nodes(channels(2061).det_node_idx);
>> det_idx = channels(2061).det_idx;
>> det_optode_idx = det_node.dets(det_idx).optode_idx;
>> det_dock_idx = enum.groups.layout.dockmap(det_node.id);
>> det_dock = enum.groups.layout.docks(det_dock_idx);
>> det_optode = det_dock.optodes(det_optode_idx).coords_3d

det_optode = 

  struct with fields:

    x: 148.5749
    y: 165.6573
    z: 48.6595
```

Et voilà.

## Data

The data structure contains the raw output data:

```
>> data

data = 

  struct with fields:

    chn_dat: [3456×978 single]
    chn_fps: 10
     chn_dt: 100
    nframes: 978
      nchns: 3456
```

The `chn_dat` field contains the raw data in an array of single precision numbers, collected at a rate of `chn_fps = 1/(chn_dt/1000)` frames per second. 

Additional fields such as channel saturation flags, accelerometry, and gyro data may be available depending upon the software and firmware versions of your system. Contact Gowerlabs to update your system.

## Layout

The layout defines the physical organisation of a LUMO group.

Let us examine the layout structure:

```
>> enum.groups.layout

ans = 

  struct with fields:

           id: '0x01338CBE'
         name: 'C005A'
    landmarks: [5×1 struct]
        docks: [1×12 struct]
      dockmap: [12×1 double]
```

The `id` is a unique identifier which matches the layout to the recorded group, though this is not enforced by software.

The docks field is an array of dock descriptors, we choose an arbitrary index:

```
>> enum.groups.layout.docks(5)

ans = 

  struct with fields:

         id: 5
    optodes: [1×7 struct]
```

The `id` field of the dock will always match with the node ID of the tile it accommodated. The optodes array contains details of each optode, as we have previously seen:

```
>> enum.groups.layout.docks(5).optodes(2)

ans = 

  struct with fields:

        name: '2'
   coords_2d: [1×1 struct]
   coords_3d: [1×1 struct]
```

Because 

The `dockmap` field is a convenience vector which maps from node/dock IDs to an index into the `docks` array. That is to say that `dockmap(4)` should return the index into the `docks` array containing the dock having an ID of 4.

Tha `landmarks`field contains a set of named positions on the group in the same 3D co-ordinate space as the optodes, this can be used for registration. The choice and naming of the landmarks is arbitrary from the perspective of this package. For example:

```
>> enum.groups.layout.landmarks(1)

ans = 

  struct with fields:

        name: 'Nasion'
   coords_3d: [1x1 struct]
```

# Bundled dependencies

### matlab-toml

lumomat bundles a [fork](https://github.com/gaetawoo/matlab-toml) of the [matlab-toml](https://github.com/g-s-k/matlab-toml) library to read and write TOML files. The package code has been re-organised to permit appropriate namespacing, and the original licence file is located at `+lumofile/+toml/LICENSE`.
