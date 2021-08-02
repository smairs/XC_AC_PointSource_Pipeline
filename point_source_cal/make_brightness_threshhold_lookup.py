import json

# Define Brightness Thresholds to use for selecting
# groups of point sources in each region that, together
# when used as calibrators, we expect to return a
# flux uncertainty of different target percentages.
# These values were derived by looking at plots of the reciprocal of the 
# ensemble signal to noise: (Flux * sqrt[N] / noise)^-1
#
#
#
#
#
#
# As of August 2, 2021 - only thresholds expect to return 5% target uncertainties are defined.

# Organise the brighness thresholds by wavelength, region, and target uncertainty in percent
brightness_threshes = {}

brightness_threshes['850'] = {}
brightness_threshes['450'] = {}

brightness_threshes['850']['IC348']   = {}
brightness_threshes['850']['NGC1333'] = {}
brightness_threshes['850']['NGC2024'] = {}
brightness_threshes['850']['NGC2071'] = {}
brightness_threshes['850']['OMC23']   = {}
brightness_threshes['850']['OPHCORE'] = {}
brightness_threshes['850']['SERPM']   = {}
brightness_threshes['850']['SERPS']   = {}
brightness_threshes['850']['DR21C']   = {}
brightness_threshes['850']['DR21N']   = {}
brightness_threshes['850']['DR21S']   = {}
brightness_threshes['850']['M17']     = {}
brightness_threshes['850']['M17SWex'] = {}
brightness_threshes['850']['S255']    = {}

brightness_threshes['450']['IC348']   = {}
brightness_threshes['450']['NGC1333'] = {}
brightness_threshes['450']['NGC2024'] = {}
brightness_threshes['450']['NGC2071'] = {}
brightness_threshes['450']['OMC23']   = {}
brightness_threshes['450']['OPHCORE'] = {}
brightness_threshes['450']['SERPM']   = {}
brightness_threshes['450']['SERPS']   = {}
brightness_threshes['450']['DR21C']   = {}
brightness_threshes['450']['DR21N']   = {}
brightness_threshes['450']['DR21S']   = {}
brightness_threshes['450']['M17']     = {}
brightness_threshes['450']['M17SWex'] = {}
brightness_threshes['450']['S255']    = {}

# For 850 microns, Low mass regions - the values are in Jy/beam to match Doug's original catalogues produced in 2017.
brightness_threshes['850']['IC348']['5.0']   = 0.3
brightness_threshes['850']['NGC1333']['5.0'] = 0.4 
brightness_threshes['850']['NGC2024']['5.0'] = 5.0
brightness_threshes['850']['NGC2071']['5.0'] = 0.5 
brightness_threshes['850']['OMC23']['5.0']   = 2.0
brightness_threshes['850']['OPHCORE']['5.0'] = 2.0 
brightness_threshes['850']['SERPM']['5.0']   = 0.8
brightness_threshes['850']['SERPS']['5.0']   = 0.9

# For 850 micron, High mass regions - the units are mJy/beam to match the catalogues Steve produced in 2021 
brightness_threshes['850']['S255']['5.0']    = 1e3
brightness_threshes['850']['M17SWex']['5.0'] = 2e3
brightness_threshes['850']['M17']['5.0']     = 1.5e3
brightness_threshes['850']['DR21S']['5.0']   = 4.5e2
brightness_threshes['850']['DR21N']['5.0']   = 4.5e2
brightness_threshes['850']['DR21C']['5.0']   = 2e3

# For 450 microns (all regions), the units are mJy/beam to match the catalogues Steve produced in 2021 
brightness_threshes['450']['IC348']['5.0']   = 900.0
brightness_threshes['450']['NGC1333']['5.0'] = 2e3
brightness_threshes['450']['NGC2024']['5.0'] = 8e3
brightness_threshes['450']['NGC2071']['5.0'] = 1.5e3
brightness_threshes['450']['OMC23']['5.0']   = 3e3
brightness_threshes['450']['OPHCORE']['5.0'] = 9e3
brightness_threshes['450']['SERPM']['5.0']   = 1.5e3
brightness_threshes['450']['SERPS']['5.0']   = 1.3e3
brightness_threshes['450']['DR21C']['5.0']   = 1.3e4
brightness_threshes['450']['DR21S']['5.0']   = 1.8e3
brightness_threshes['450']['DR21N']['5.0']   = 2e3
brightness_threshes['450']['M17']['5.0']     = 2e4
brightness_threshes['450']['M17SWex']['5.0'] = 1e4
brightness_threshes['450']['S255']['5.0']    = 4e3

with open('point_source_cal/bright_threshes.json','w') as outfile:
    json.dump(brightness_threshes,outfile) 
