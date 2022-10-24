# lumomat

## v1.4.1
  - do not use the memory command on platforms upon which it is not supported
  - update documentation on MPU data output ordering, extended metadata

## v1.4.0
  - fix bug reading certain LUFR files with embedded layouts
  - add bad frame markers to data and output files

## v1.3.0
 - read LUFR files v3
 - support embeeded layouts in LUFR files
 - support hyperscanning in LUFR files

## v1.2.0

- remove 'mne-nirs' style output, default now supported upstream
- add string to character array input validation

## v1.1.0

 - read LUMO files v0.5
 - add `merge_layout` function to integrate layouts into LUMO files
 - saturation data is now chunked and compressed
 - SNIRF extended metadata now stored in root of the file, allowing validation

## v1.0.0

Initial release
 - read LUMO files v0.1 -> v0.4
 - read LUFR files v1 -> v2
 - write NIRS files (including DOT-HUB style)
 - write SNIRF files, compliant with exceptions confirmed by spec. authors
 - basic data exploration




