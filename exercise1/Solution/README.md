# High Performance Computing final project - Exercise 1

## Student's Info

| Name | Surname | Student ID | UniTs email | Personal email | Master course |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Piero | Zappi | SM3600004 | `piero.zappi@studenti.units.it` | `piero.z.2001@gmail.com` | SDIC |

## Folder Organization

This folder is organized with the following structure:

```bash
.
├── Compile_OSU
│   ├── compilation.err
│   ├── compilation.out
│   ├── compilation_EPYC.sh # Bash script to compile the OSU benchmark on the EPYC nodes of the ORFEO cluster
│   └── osu-micro-benchmarks-7.3 # OSU benchmark
├── Report_and_slides # Final report and slides for the presentation
│   ├── slides_ex1
│   └── Zappi_ex1_report.pdf
├── README.MD # This file
├── barrier # Bash scripts to automate the data collection process and notebook for the analysis of the results
│   ├── Results # Data collected on the ORFEO cluster (csv files)
│   │   ├── barrier_bruck.csv
│   │   ├── barrier_default.csv
│   │   ├── barrier_doublering.csv
│   │   ├── barrier_linear.csv
│   │   └── barrier_tree.csv
│   ├── barrier_bruck.sh
│   ├── barrier_data_analysis.ipynb
│   ├── barrier_default.sh
│   ├── barrier_doublering.sh
│   ├── barrier_linear.sh
│   └── barrier_tree.sh
└── bcast # Bash scripts to automate the data collection process and notebook for the analysis of the results
    ├── Results # Data collected on the ORFEO cluster (csv files)
    │   ├── bcast_basiclinear.csv
    │   ├── bcast_binarytree.csv
    │   ├── bcast_chain.csv
    │   └── bcast_default.csv
    ├── bcast_basiclinear.sh
    ├── bcast_binarytree.sh
    ├── bcast_chain.sh
    ├── bcast_data_analysis.ipynb
    └── bcast_default.shcast

```

## Contents

This folder contains the solution for exercise 1, which requested to compare and assess different openMPI algorithms for collective operations and to develop models to predict their latencies; specifically, the broadcast and the barrier operation were chosen.
Further details can be found in the **Zappi_ex1_report.pdf** file.
