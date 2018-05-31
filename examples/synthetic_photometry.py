# import the HCSS point source simulation class
from simulation.photometric.point_source import SimPointSource
# import the simulation plot environment
from out.plot.plot import Plot

# set the path to the synthetic photometry and to the save place
mist_path = './mist/'
save_path = './save/'

# Create a point source simulation setup with the default settings
# take the synthetic photometry from the mist_path place
sps = SimPointSource(mist_path)

# save the simulation parameters in this directory
sps.parameters.save_html('./')

# start the simulation and derive the photometry for 1000 clusters
sps.start(1000)

# save the results
sps.save(save_path)

# create a plot environment for the simulation results
p = Plot(sps)

# iterate over some time_id's
# a time id is just the position of the time in the time-array
for time_id in [5, 10, 20]:

    # plot the magnitude distribution of Pan-STARRS r and 2MASS Ks
    p.brightness_distribution('PS_r', time_id)
    p.brightness_distribution('2MASS_Ks', time_id)

    # plot the color distribution (first band minus second band)
    p.color_distribution('PS_g', 'PS_r', time_id)

    # plot a color-color distribution (first band minus second band and
    # third band minus fourth band)
    p.color_map('PS_g', 'PS_r', 'PS_r', 'PS_i', time_id)

    # plot a color magnitude distribution (first band - second band and the second
    # band will be used for the magnitude axis)
    p.color_mag_map('PS_g', 'PS_r', time_id)
    p.color_mag_map('WISE_W1', '2MASS_J', time_id)

    # plot a color magnitude distribution (first band - second band,
    # instead of the default second band as the magnitude axis
    # use the band3 value as the magnitude axis
    p.color_mag_map('Bessell_V', 'Bessell_I', time_id, band3='Bessell_V')
    p.color_mag_map('Bessell_V', '2MASS_Ks', time_id, band3='Bessell_V')
