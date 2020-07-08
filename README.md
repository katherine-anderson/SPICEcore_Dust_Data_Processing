# SPICEcore Dust Data Processing
**This repository contains code for processing of the 54.3 ka South Pole Ice Core (SPICEcore) dust flux and dust size distribution records.**

Code was developed in Python version 3.8. Go [here](https://docs.anaconda.com/anaconda/navigator/install/) to install the Anaconda Navigator. Code was originally developed with the [Spyder](https://www.spyder-ide.org/) application. Code and data must be downloaded before running. We recommend saving all code files to one folder and all data files to another folder.

## Attribution
Data processing code was developed in close collaboration with Aaron Chesler for the following dissertation:

Anderson, Katherine L. (2020). *Atmospheric Dynamics during the Abrupt Climate Change Events of the Last Glacial Period.* Location: Dartmouth College.

## Code Files
- "Complete_SPICEcore_Dust_Processing_Final.py"
  - Master data processing script. All code can be run from this file.
  - Dust processing has 2 phases
    - Removal of melting errors 
    - Removal of contamination
  - To run both phases of data processing, enter 'Y' on the 1st prompt
  - Select 'N' on the 1st prompt if you ran the Phase 1 script already and have a "CLEANED_CFA_Phase1..." dataset
  
- Phase 1
  - Occurs in "SPICEcore_Dust_Phase1_Processing_Final.py"
  - Counts and removes melting errors
  - Applies the SP19 timescale (Winski et al., 2019)
  - Adds descriptive columns
  - Calculates dust metrics
  - Saves cleaned data ("Cleaned_CFA_Phase1...")
  
- Phase 2
  - Occurs in "Complete_SPICEcore_Dust_Processing_Final.py"
  - Preserves data during volcanic events and dust events
  - Removes outliers
  - Removes remaining manually-identified issues
  - Prints summary statistics
  - Saves removed data ("Bad_CFA...") and cleaned data ("Cleaned_CFA_Phase2...")

- "SPICEcore_Dust_Processing_Functions_Final.py"
  - Contains functions used in Phase 1 and Phase 2 data processing
  
## Data Files
Data are stored in a separate repository and are accessible at [this link](https://rcweb.dartmouth.edu/homes/f003qyw/).
