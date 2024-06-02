# radtk : A toolkit for working with RAD files

### About 

`radtk` is collection of various tools for working with RAD files.  `radtk` is very much a work in progress, so suggestions for new tools, features and capabilities is very much welcome, 
as are contributions and PRs.

### Basic usage

`radtk` is a program, written in [`rust`](https://www.rust-lang.org/), and making use of [`libradicl`](https://github.com/COMBINE-lab/libradicl) for performing various operations on RAD files. 
It consists of a number of different sub-commands.  The current sub-commands are listed below:

 - `cat` : Concatenates a series of input RAD files into an output RAD file.  The input files must have compatible headers (i.e. be built on the same referece, and contain the 
 same tag set).

 - `view` : Output a RAD file into a textual (currently JSON) format.  This can be useful for inspecting the contents of a RAD file or specifically examining or searching for some record 
 or quantity of interest.
