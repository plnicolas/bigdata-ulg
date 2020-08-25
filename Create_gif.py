import imageio
import os
images = []

path = '/CECI/home/ulg/mast/eivanov/Big_Data/Pics_for_animation/'
d = sorted(os.listdir(path))
for filename in d:
	print (filename)
	images.append(imageio.imread(path+filename))
kargs = { 'duration': 1 }
imageio.mimsave('Results/Relative_trend.gif', images, 'GIF', **kargs)
