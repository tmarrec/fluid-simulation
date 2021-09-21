import os
import glob
from xml.dom import minidom

ssp = 2048
width = 1024
height = 1024
plydir = 'save-ply'
resdir = 'res'

os.system('mkdir -p '+resdir)
meshesdir = sorted(glob.glob(plydir+'/*.ply'), key=os.path.getmtime)
for mesh in meshesdir:
    xml = "<?xml version='1.0' encoding='utf-8'?>\
    <scene version=\"0.6.0\">\
    <integrator type=\"volpath_simple\">\
            <integer name=\"maxDepth\" value=\"24\"/>\
    </integrator>\
    <sensor type=\"perspective\">\
            <float name=\"focusDistance\" value=\"6\"/>\
            <float name=\"fov\" value=\"60\"/>\
            <string name=\"fovAxis\" value=\"x\"/>\
            <transform name=\"toWorld\">\
                    <lookAt target=\"32.0, 16.0, 32.0\" origin=\"120.0, 90.0, 45.0\" up=\"0, 1, 0\"/>\
            </transform>\
            <sampler type=\"stratified\">\
                    <integer name=\"dimension\" value=\"8\"/>\
                    <integer name=\"sampleCount\" value=\""+str(ssp)+"\"/>\
            </sampler>\
            <film type=\"hdrfilm\">\
                    <boolean name=\"banner\" value=\"false\"/>\
                    <integer name=\"width\" value=\""+str(width)+"\"/>\
                    <integer name=\"height\" value=\""+str(height)+"\"/>\
                    <string name=\"pixelFormat\" value=\"rgb\"/>\
                    <rfilter type=\"gaussian\"/>\
            </film>\
    </sensor>\
    <emitter type=\"sky\">\
            <float name=\"scale\" value=\"12\"/>\
            <float name=\"turbidity\" value=\"5\"/>\
    </emitter>\
    <shape type=\"rectangle\">\
            <transform name=\"toWorld\">\
                    <rotate x=\"1\" angle=\"-90\"/>\
                    <scale value=\"1000\"/>\
                    <translate y=\"0\" />\
            </transform>\
            <bsdf type=\"diffuse\">\
                    <texture name=\"reflectance\" type=\"checkerboard\">\
                            <float name=\"uscale\" value=\"50\"/>\
                            <float name=\"vscale\" value=\"50\"/>\
                    </texture>\
            </bsdf>\
    </shape>\
    <shape type=\"ply\">\
            <string name=\"filename\" value=\""+mesh+"\"/>\
            <bsdf type=\"dielectric\">\
                    <string name=\"intIOR\" value=\"water\"/>\
                    <string name=\"extIOR\" value=\"air\"/>\
            </bsdf>\
    </shape>\
    </scene>\
    "
    n = mesh.split('/')[-1].split('.ply')[0]
    file = open(str(n)+'.xml', 'w')
    file.write(minidom.parseString(xml).toprettyxml(indent="    "))
