
import matplotlib.pyplot as plt
import numpy as np
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import os

x = np.linspace(-2, 2, 200)

duration = 2

savedir = '/home/laojin/my_lat/moviepy/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

fig, ax = plt.subplots()
def make_frame(t):
	print(t)
	ax.clear()
	ax.plot(x, np.sinc(x**2) + np.sin(x + 2*np.pi/duration * t), lw=3)
	ax.set_ylim(-1.5, 2.5)
	return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration=duration)
animation.write_videofile(savedir+"matplotlib.mp4", fps=24,codec = 'h264')
animation.write_gif(savedir+"matplotlib.gif", fps=24)
animation.close()






