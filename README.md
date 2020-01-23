# FalkorFactor
A repository for code used to process and analyze data from the Falkor cruise. 


## Cruise overview
Metadata and some nice visualizations of the cruise are available in the MetaDataVis folder.

## Data formats
Raw files are stored on the Ingalls QE drive, and converted to .mzML 
via the RunMsconvert batch file (but ignored by Git).

MSDIAL requires .abf files, so those are converted from .mzML and stored in the MSDIAL folder along with all 
the other MSDIAL side effects (but ignored by Git).

XCMS can operate directly on the .mzML files, so the XCMS folder is largely focused on scripts and processing information as well as
intermediate objects.