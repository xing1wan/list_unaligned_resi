Author: Xing Wan Date: 2024-03-15

Inspired by rmsdByRes Zhenting Gao on 7/28/2016 and rmsdCA Yufeng Tong 2020-01-29

# A PyMOL script to calculate alpha carbon distances between two aligned protein structures

This is a `PyMOL` script to calculate the Euclidean distance between of two aligned alpha carbons (CAs) and to highlight unaligned residues.

The script will select the CAs from the reference and target protein stuctures, and align them at a default cutoff 2Å. The aligned residues will be coloured grey, and unaligned residues from the reference and target will be red and blue, respectively.

## Installation

Save the script to any folder accessible in `PyMOL`.

## Usage
This script has been tested on `PyMOL` version 2.7.

1. Load two structures in the `PyMOL` program.
2. `Set_name` of your two stuctures to more convenient ones.
3. `File` → `Run` this script.
4. Type in the command line, e.g., `list_unaligned_resi('/ref_structure//A', '/tgt_structure//A', 2.0)` or `list_unaligned_resi('ref', 'tgt', 2)`
5. There is a choice of show the aligned and unaligned residues in surface, `uncomment` lines 56, 60, 65, 69 to enable `surface`.

## Possible updates
1. To allow cutoff values set by the user
2. To add other visualisation choices
