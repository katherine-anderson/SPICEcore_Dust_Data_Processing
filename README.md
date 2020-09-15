# SPICEcore Dust Data Processing
**This repository contains code for processing of the 54.3 ka South Pole Ice Core (SPICEcore) dust flux and dust size distribution records.**

Code was developed in Python version 3.8. Go [here](https://docs.anaconda.com/anaconda/navigator/install/) to install the Anaconda Navigator. Code was originally developed with the [Spyder](https://www.spyder-ide.org/) application. Code and data must be downloaded before running. We recommend saving all code files to one folder and all data files to another.

## Attribution
Data processing code was developed in close collaboration with Aaron Chesler for the following dissertation:

Anderson, Katherine L. (2020). Atmospheric Dynamics during the Abrupt Climate Change Events of the Last Glacial Period, (master's thesis). Retrieved from ProQuest Dissertations Publishing. (Order No. 28089979). Hanover, NH: Dartmouth College.
Available at [this link] (https://drive.google.com/file/d/1THwkHyc8u7LrJVMey_zb3pkgKeHBiS8l/view?usp=sharing).

## Data Files
Input and output data are stored in a separate repository and are accessible at [this link](https://rcweb.dartmouth.edu/homes/f003qyw/).

## Code Files
- *"Complete_SPICEcore_Dust_Processing.py"*
  - **Master data processing script. All data processing can be run from this file.**
  - Dust processing has 2 phases
    - Removal of melting errors 
    - Removal of contamination
  - To run both phases of data processing, enter 'Y' on the 1st prompt
  - Select 'N' on the 1st prompt if you ran the Phase 1 script already and have a "Cleaned_CFA_Phase1..." dataset
  - Code will ask the user for paths to the code and data folders
  
- Phase 1 data cleaning
  - Occurs in *"SPICEcore_Dust_Phase1_Processing.py"*
  - Code will ask the user for paths to the code and data folders
  - Counts and removes melting errors
  - Applies the SP19 timescale (Winski et al., 2019)
  - Adds descriptive columns (see below)
  - Calculates dust metrics (see below)
  - Saves cleaned data (*"Cleaned_CFA_Phase1..."*)
  
- Phase 2 data cleaning
  - Occurs in *"Complete_SPICEcore_Dust_Processing.py"*
  - Preserves data during dust events
  - Gives user the option to preserve data during volcanic events
  - Removes outliers
  - Removes remaining manually-identified issues
  - Prints summary statistics
  - Saves removed data (*"Bad_CFA..."*) and cleaned data (*"Cleaned_CFA_Phase2..."*)

- Functions used in Phase 1 and Phase 2 data cleaning
  - In *"SPICEcore_Dust_Processing_Functions.py*"
  - Other script files automatically read function definitions from here
  
- Columns added to the raw data during data cleaning
  - "AgeBP": Age (years before 1950) based on SP19 timescale (Winski et al., 2019)
  - "Break?": "True" if data fall within specified depth range of core breaks
  - "New Break?": "True" for the first row occurring within each core break window (can be used to count core breaks, for example)
  - "Volcanic Event?": "True" if data fall within specified time range of volcanic events
  - "New Volcanic Event?": "True" for the first row occurring within each volcanic event window (can be used to count volcanic events, for example)
  - "Dust Event?": "True" if data fall within designated depth ranges of dust events which were visible during melting
  - "Sum 1.1-12": Particle number concentration (# of particles ≥1.1 µm diameter/µL)
  - "CPP": Coarse particle percentage (particles ≥4.5 µm / particles ≥1 µm * 100; after Koffman et al., 2014)

- "Old Scripts" folder: script archive, not required for data processing
  

