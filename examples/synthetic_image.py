# import the HCSS point source simulation class
from simulation.photometric.point_source import SimPointSource
# import the HCSS image simulation class
from simulation.image.image import ImageSim
# import the image SExtractor interaction class
from simulation.image.img_processing import ImgSimSExtractor

mist_path = './mist/'
# create a synthetic photometry simulation object
sps = SimPointSource(mist_path)
# select a metallicity of 0.0
sps.directory = '0_0'
# create an image simulation object with a radius of 0.1 pc and a distance module of 15
img = ImageSim(sps, 0.1, distance_module=15)

# calculate the 2d projected position
img.calc_pos()
# create and store the SkyMaker star input list with the first survey
img.create_skymaker_output('./', sps.parameters.surveys[0], 'PS_r')

# load the object for the SExtractor data
img_sext = ImgSimSExtractor(img)

# plot the simulated magnitude vs the extractor magnitude from the simulated image
img_sext.plot_magnitudes()
