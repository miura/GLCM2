from emblcmci.glcm import GLCMtexture
from ij import IJ
imp = IJ.openImage('/Volumes/D/Julia20130201-/NucleusSegmentationStudy/20130408/grayscale12.tif')

glmc = GLCMtexture(1, 45, True, False)
ip = imp.getStack().getProcessor(1)
ip8 = ip.convertToByte(False)
glmc.calcGLCM(ip8)
resmap = glmc.getResultsArray()
pa = glmc.paramA
for p in pa:
    print p, resmap.get(p)
