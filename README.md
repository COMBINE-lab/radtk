# radtk : A toolkit for working with RAD files

A collection of various tools for working with RAD files.  `Radtk` is very much a work in progress, so suggestions for new tools, features and capabilities is very much welcome, 
as are contributions and PRs.

### Current tools

 - `cat` : Concatenates a series of input RAD files into an output RAD file.  The input files must have compatible headers (i.e. be built on the same referece, and contain the 
 same tag set).
 - `viea` : Output a RAD file into a textual (currently JSON) format.  This can be useful for inspecting the contents of a RAD file or specifically examining or searching for some record 
 or quantity of interest.
