# Some functions related to PyMol

## Background

While working to both analyze and visualize some similar structures, I decided to converte some code into functions that can be reused. No guarantees that they will work for your purposes.

I've tried to document the logic with in-line comments which hopefully will provide you with guidance when you need to make modifications and improve on what I've done.

## Current functions

- `draw_caver_mesh`: This function attempts to maximize your Caver cavity visualization by finding spheres from other clusters that are near a main cluster you've identified manually. It then draws a mesh representation for this final object.
- `get_avg_sc_b`: Gets average sidechain B-factor given a PyMol selection as a string.
- `bfac_comp`: Relies on `get_avg_sc_b` to get the average sidechain B-factor across multiple identical models and calculates the mean and standard deviation (and normalized forms) across those models at each position. The normalization, for whatever reason, is across the residues you picked and not across the entire structure. I need to fix that.

## Conda environment .yml

I've included the yml file from `conda env export` for reference. There are more packages there than you need. Currently you only need Jupyter, pandas, and opensource PyMol to run the scripts I have included.
