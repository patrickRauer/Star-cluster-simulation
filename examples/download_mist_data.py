from mist.mist_download import download_all

# set the boundaries for the download

# list with the required metallicities in Fe/H
metallicities = [-0.2, 0.0, 0.2, 0.4, 0.5]
# set the minimal mass
min_mass = 0.1
# set the maximal mass
max_mass = 2.0
# set the mass resolution
mass_step = 0.001
# set the directory where the downloaded data will be stored
directory = './mist/0_2/'

# start the download program
download_all(min_mass, max_mass, mass_step, metallicities,
             directory)
