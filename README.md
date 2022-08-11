# vertical_growth_biofilms

This is all the code utilized in the analysis for "Vertical
growth dynamics of biofilms". 

- First we build a `.csv` database for all the timelapse data using `database_creation.jl`, many experimental metrics are calculated here.
- Profiles are then fitted to multiple models in `/fitting/allstrainfit.jl`
- Code for each figure is in `/figs/scripts/`.

Original `.datx` files from the Zygo Interferometer are available upon request.