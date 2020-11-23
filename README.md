These are analysis scripts written for analyzing IFU datacubes obtained through FIFI spectroscopy done by SOFIA, NASA's space telescope. The scripts perform the following functions -

1. Fluxmaps.ipynb - this creates a spectral density heatmap across the observed spaxels of the galaxy for all available wavelengths

2. Wavelength_Profile_Central.ipynb - this creates a spectral density plot against wavelength for a central n_x x n_y spaxel^2 window (user specified) of the galaxy

3. x_[y].ipynb - this is an analysis pipeline that computes the flux value of line "y" in its observed spectrum from galaxy "x"