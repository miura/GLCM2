importClass(Packages.emblcmci.glmc.GLMCtexture);

for (var i = 1; i < 20; i++){
	g = new GLMCtexture(IJ.getImage(), i, 45, true, false);
	g.calcGLMC();
	ht = g.getResultsArray();
	IJ.log(ht.get("Contrast")) ;
	g.writetoResultsTable();
}