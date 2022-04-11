# lumomat
MATLAB tools for LUMO data

# Installation

1. Download the latest version from the [release page](https://github.com/Gowerlabs/lumomat/releases)
2. Unzip in your favourite directory
3. Add the  `lumomat` directory to your MATLAB path

# Quickstart

Load a `.lumo` file:

```
data = LumoData(filename)
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
       src_coord_2d: [94.3617 75.6020]
       src_coord_3d: [87.3064 179.5761 110.9974]
        det_node_id: 2
            det_idx: 2
    det_optode_name: '2'
       det_coord_2d: [102.0520 81.9720]
       det_coord_3d: [79.4308 175.0386 115.1332]
```

Find all of the channel indices where the source node ID is 1, the detecot node ID is 8, and the wavelength is 850nm:

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

And the channel data for this channel

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

Convert a `.lumo` file to a SNIRF file:

```
data.write_NIRS(filename);
```

Conver a `.lumo` file to a NIRS file:

```
data.write_SNIRF(filename);
```



# Overview

LUMO is a high-density, wearable, and modular instrument for diffuse optical tomography (DOT) and (functional) near-infrared spectroscopy ((f)NIRS). 

In principle, the data produced in a DOT experiment is striaghtforward: a single mesaurement of intensity is made across every channel of the system many times per second, resulting in a big matrix of data. In practice, there is a signficant amount of book-keeping required in order to make sense of the relationship between that big matrix of numbers and the physical geometry of the system.

This is particularly true for LUMO, owing to its modular nature. The purpose of this package is to allow one to access, interpret, and export the data produced by LUMO, for the purposes of further analysis.

# Features

 - High-level stable object based API for interacting with lumo output 
 - Low-level functional API for custom use-cases
 - Load .lumo files into a standardised format for further analysis
 - Write data to NIRS or SNIRF format


# Nomenclature

Various terms are used throughout this guide and in the package itself. Many will be familiar to practitioners in (f)NIRS and DOT, others are specific to LUMO, and the most dangerous are ambiguous.

 - *Tile* / *Node* / *Module*: a LUMO tile is a single optoelectronic module containing a set of optical sources, detectors, motion processing units, and so on. Lumo tiles are plugged into docks to form a DOT system. In the logical representation of the system, a tile is called a 'node' (all tiles are nods, but not all nodes are tiles, but for most practical situations they are synonymous). A node is the same concept as a 'module' in the parlance of the SNIRF file format.
 - *Dock*: a dock is the recepticle into which a LUMO tile is installed. Multiple docks are normally connected together in a piece of headgear, or a patch, which is applied to the subject of the experiment.
 - *Group* / *Cap*: most fNIRS and DOT experiments involve mesaurements of the heamodynamics of the brain, and thus most systems are formed of headgear called caps into which multiple docks are installed. But there are many other interesting targets which can be measured, so the generic name for a gorup of docks which form a single mesaurement system is, imaginatively, a 'group'.
 - *Optode*: an optode is a single location from which light is emitted by one or more sources, and/or detected by one or more detectors. Conceptually there is a one-to-one mapping between the optode on a tile, and the optode on a dock: this is the link between logical sources and physical locations.
 - *Source*: LUMO considers a source to mean a single optical emitter at a specific wavelength. There are usually multiple sources at different wavelengths co-located at a single optode in order to enable spectroscopy. *This meaning may differ from some file formats and analysis packages which think of a source as an optode, with multiple wavelengths.*
 - *Detector*: a detector is a device which measures optical intensity. Like sources, detectors are associated with a particular optode.
 - *Channel*: a channel is formed from a single source-detector pair.
 - *Layout*: a layout describes the docks present in a group, and each of the optodes which belong to the dock. When a cap is built, a template layout is created which provides the positions of all docks in a group as measured on a suitable phantom. Subject specific layouts might be mesaured by the user (and may consist of a subset of the docks in the group).
 - *Frame*: a frame is a measurement of all active channels made at a single point in time.


# Using matlumo







## Notes

 - If your .lumo file does not contain a template layout file a warning will be issued and it will not be possible to write the data to different output formats. See the 'Layouts' section of this document to resolve this problem.

!!!
-- Add layout details to this section ---
!!!


## Local and global indexing

The canonical enumeration of a LUMO system uses local indexing. This means that channels are defined by a 4-tuple: 

*ch = (source node, source index, detector node, detector index)*

This representation permits complete flexibility in the description of channels in the system, and maintains the link between channels, nodes, and docks, which is required in a modular system.

Many data formats and analysis tools for (f)NIRS and DOT use an alternative format to index the channels of the system. The most common choice is the global spectroscopic scheme, in which a channel is defined by a 3-tuple:

*ch = (source optode, detector optode, source wavelength)*

This indexing scheme is less flexible, but it is sane because:

 - most experiments involve spectoscopy, so it is a fair assumption that all source optodes will have associated with them sources at every wavelength
 - most systems do not collocate sources and detectors in the same optode, so it is reasonable to assign an optode to be either belonging to sources, or detectors

The canonical enumeration of a LUMO system can usually be converted to global indexing, allowing the data to be exported to a number of formats.

## NIRS output

Additional fields

## SNIRF output

Additional fields




# Low-level functional API

In most circumstances, users should choose the high-level object based API in order to access their data. The high-level wrapper is implemented on-top of a low-level functional API, namespaced in MATLAB packages. 

No gaurantees are made regarding the stability of the low-level API, so users must take account of versioning when calling the low-level API.





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

Whilst loading the file, information will be printed regarding the contents of the file. Some `.lumo` files might not contain an embedded layout file. If this is the case a warning will be printed:

```
Warning: The sepcified LUMO file does not contain an embedded layout file, and no layout has been specified when calling this function. The returned enumeration will lack layout information, and it will not be possible to convert this file to formats which require a layout. Specify an appropriate layout file to supress this warning, or copy an appropriate layout to the .LUMO folder in order for it to be used as an automatic fallback. 
```

Follow the guidance, or read more about Layouts in a later section.

## Enumeration, Enumeration, Enumeration

When you load a `.lumo` file, `lumomat` constructs a canonical enumeration which provides a full description of the system. Navigating the enumeration allows the data to be interpreted and manipulated into the form required for further analysis.

The file we have loaded was recorded on a system using 12 tiles, resulting in 3456 channels. The enumeration contains all the information required to interpret the raw data. We can examine its contents:

```
>> enum

enum = 

  struct with fields:

       hub: [1×1 struct]
    groups: [1×1 struct]
```

The `hub` field contains information about the system that is not typically used during data analysis. The `groups` field is an array of structures, one for each of the groups connected to the LUMO hub. Since `.lumo` files only record a single group at a time, it is a scalar structure and no indexing is required in this case:

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

The `nodes` array contains information about each node that formed the system, incuding their sources and detectors. For example, the eigth entry:

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

This node happens to have an `id` of 8, as well as being located in the eigth element of the array. This is because the system was fully popoulated. In general, the node ID and the index in the array do not correspond. We can see from the dimensionality of the `srcs` field that this node has six sources. Let's examine source one:

```
>> src_idx = 1;
>> enum.groups.nodes(node_idx).srcs(src_idx)

ans = 

  struct with fields:

            wl: 735
    optode_idx: 5
```

And we can extract information about all of them, for example, the unqiue wavelengths:

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

The name of the optode matches that used in the user manual for the LUMO system, it is only used for display purposes. This optode is of type 'S' which means that it is only associated with optical sources. This is of signficance when we discuss global indexing.

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
    coord_2d: [1×1 struct]
    coord_3d: [1×1 struct]
```

And finally determine its location:

```
>> optode.coord_3d

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

Since we have already developed significant affection for the first source on the eigth node, let's pick a channel for which it provides the illumination:

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

The source node and source indices should not come as a suprise. The wavelength of this channel is:

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
>> det_optode = det_dock.optodes(det_optode_idx).coord_3d

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

Whilst the majority of the data in the enumeration comes from the system itself, the layout information is copied from a file which is built when the particular group is manufactured, based on measurements made on a suitable phantom. We refer to this layout as the *default* layout.

*Note: some versions of the LUMOview software do not embed a layout file in the `.lumo` output file. If no layout file is embedded, a warning will be printed and the resultant layout structure will be empty. In this instance the `read_lumo` function will fall back to load any file named `layout.json` within the `.lumo` directory, so it is possible to copy, rename, and move the Gowerlabs' supplied layout file to restore functionality. Alternatively, the layout file may be specified as described for the provision of subject-specific layouts.*

Let us examine the layout structure:

```
>> enum.groups.layout

ans = 

  struct with fields:

          uid: 20155582
    landmarks: [5×1 struct]
        docks: [1×12 struct]
      dockmap: [12×1 double]
```

The `uid` is a unique identifier which matches the layout to the recorded group, though this is not enforced by software.

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
    coord_2d: [1×1 struct]
    coord_3d: [1×1 struct]
```

Because 

The `dockmap` field is a convenience vector which maps from node/dock IDs to an index into the `docks` array. That is to say that `dockmap(4)` should return the index into the `docks` array containing the dock having an ID of 4.

Tha `landmarks`field contains a set of named positions on the group in the same 3D co-ordinate space as the optodes, this can be used for registration. The choice and naming of the landmarks is arbitrary from the perspective of this package. For example:

```
>> enum.groups.layout.landmarks(1)

ans = 

  struct with fields:

    name: 'Nasion'
       x: 85.5734
       y: 203.8798
       z: 17.4744
```


It is often the case that for the purposes of data analysis, practitioners will wish to provide subject-specific layout information measured during an experiment. This information can be used instead of the default layout:

```
>> [enum, data, events] = lumofile.read_lumo('sample.lumo', 'layout', custom_layout);
```

where `custom_layout` can be one of:
 - a structure matching that described here, or, 
 - an JSON file following the syntax of the default layout file (contact Gowerlabs for details)

For a layout file to be usable, it must contain the loctions of all docks for which nodes were present. It is acceptable therefore that only a subset of the docks are recorded, corespondig to those docks where nodes are present.



# Bundled dependencies

### matlab-toml

lumomat bundles a [fork](https://github.com/gaetawoo/matlab-toml) of the [matlab-toml](https://github.com/g-s-k/matlab-toml) library to read and write TOML files. The package code has been re-organised to permit appropriate namespacing, and the original licence file is located at `+lumofile/+toml/LICENSE`.