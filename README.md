# BTR (Ball and Turn Randomization) R Project

This R project contains code for implementing and analyzing Ball and Turn randomization procedures for clinical trials.

## Project Structure

```
BTR/
├── BTR.Rproj           # R Project file
├── README.md           # This file
├── code/               # R source code
│   ├── constructBT_v1.R    # Main BT construction functions
│   ├── exeWBT.R            # Execute Weighted BT simulations
│   ├── processWBT.R        # Process WBT results
│   ├── plotBT_work.R       # Plotting functions
│   └── try_constructWBT_v1.R # Experimental code
├── WBTsimRes/          # Simulation results output
└── code.zip            # Archived code
```

## Getting Started

1. Open `BTR.Rproj` in RStudio to activate the project environment
2. Install required packages: `install.packages(c("Rglpk", "here"))`
3. Source the main functions: `source(here("code", "constructBT_v1.R"))`

## Key Files

- `constructBT_v1.R`: Contains the main functions for constructing Ball and Turn randomization designs
- `exeWBT.R`: Executes Weighted Ball and Turn simulations
- `processWBT.R`: Processes and analyzes WBT simulation results

## Dependencies

- `Rglpk`: Linear programming solver
- `here`: Path management for reproducible file paths
