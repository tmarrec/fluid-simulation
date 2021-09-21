import os
import glob
import numpy as np
import mitsuba
mitsuba.set_variant('packet_rgb')
from mitsuba.core import Thread
from mitsuba.core.xml import load_file
from mitsuba.core import Bitmap, Struct

filename = 'scene-liquid.xml'
#meshnames = 'save-ply/*.ply'
meshnames = '../build/result/*.ply'
tempdir = 'temp'
resdir = 'images'
res = 1024
spp = 64
cubescale = 1024

eyestart = np.array([160.0, 160.0, 80.0])
eyeend = np.array([10, 0, 0])

Thread.thread().file_resolver().append(os.path.dirname(filename))

def render():
    print('Rendering')
    os.system('mkdir -p '+tempdir)
    meshesdir = sorted(glob.glob(meshnames), key=os.path.getmtime)
    iteration = 0
    eyediff = (eyeend - eyestart)/len(meshesdir)
    eyepos = eyestart
    for mesh in meshesdir:
        #mesh = '../build/result/156.ply'
        eyepos += eyediff
        eyepos = np.array([230.0, 300.0, -1.0])
        xmlpos = str(eyepos[0])+','+str(eyepos[1])+','+str(eyepos[2])
        scene = load_file(filename, meshname=mesh, eyepos=xmlpos, spp=spp, res=res, cubescale=cubescale)

        sensor = scene.sensors()[0]
        scene.integrator().render(scene, sensor)

        film = sensor.film()

        img = film.bitmap(raw=True).convert(Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True)
        img.write(tempdir+'/'+str(iteration)+'.png')
        #film.set_destination_file(tempdir+'/'+str(iteration)+'.exr')
        #film.develop()
        iteration += 1
        #exit()


def video():
    print('Video rendering..')
    print(os.popen('ffmpeg -y -r 60 -f image2 -s 1024x1024 -i '+tempdir+'/%d.png -vcodec libx264 -crf 8 -pix_fmt yuv420p render.mp4').read())
    os.system('cvlc --loop render.mp4')

render()
#video()

