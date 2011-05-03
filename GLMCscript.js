importClass(Packages.emblcmci.glmc.GLMCtexture);
g = new GLMCtexture(IJ.getImage(), 2, 45, true, true);
glmc = g.calcGLMC();
ht = g.getResultsArray(glmc);
IJ.log(ht.get("Contrast")) ;