import os
import glob
import numpy as np
import mitsuba
mitsuba.set_variant('packet_rgb')
from mitsuba.core import Thread
from mitsuba.core.xml import load_file
from mitsuba.core import Bitmap, Struct

filename = 'scene.xml'
volnames = '../build/result/*.vol'
tempdir = 'temp'
resdir = 'images'
res = 512
spp = 2048
cubescale = 4

eyestart = np.array([14.0, 2.5, 5.5])
eyeend = np.array([-3.5, -2.75, 10])

Thread.thread().file_resolver().append(os.path.dirname(filename))

def render():
    print('Rendering')
    os.system('mkdir -p '+tempdir)
    voldirs = glob.glob(volnames)
    iteration = 0
    eyediff = (eyeend - eyestart)/len(voldirs)
    eyepos = eyestart
    for vol in voldirs:
        eyepos += eyediff
        xmlpos = str(eyepos[0])+','+str(eyepos[1])+','+str(eyepos[2])
        scene = load_file(filename, volname=vol, eyepos=xmlpos, spp=spp, res=res, cubescale=cubescale)

        sensor = scene.sensors()[0]
        scene.integrator().render(scene, sensor)

        film = sensor.film()

        img = film.bitmap(raw=True).convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True)
        img.write(tempdir+'/'+str(iteration)+'.png')
        #film.set_destination_file(tempdir+'/'+str(iteration)+'.exr')
        #film.develop()
        iteration += 1


def video():
    print('Video rendering..')
    print(os.popen('ffmpeg -y -r 60 -f image2 -s 800x800 -i '+tempdir+'/%d.png -vcodec libx264 -crf 8 -pix_fmt yuv420p render.mp4').read())
    os.system('cvlc --loop render.mp4')

render()
video()

