/** GCLM
Script for analyzing stack and printing out results (faster than running the plugin interatively)
Kota Miura (miura@embl.de) 
*/
importClass(Packages.emblcmci.glcm.GLCMtexture);
importClass(Packages.java.util.HashMap);

imp = IJ.getImage();
if (imp.getNSlices() == 1) {
	IJ.log("Need a Stack")
} else {
	gl = GLCMtexture();
	gl.d = 30;	//gap distance
	gl.phi = 0;
	gl.rt_reset = false;
	var slices = imp.getNSlices();
	var glcmA = [];
	var pA = [];
	var paranamesA = gl.paramA;
	for (var i = 0; i< slices; i++){
		var ip = imp.getStack().getProcessor(i+1);
		glcmA[i] = gl.calcGLMC(ip); // returned is double[][]
		pA[i] = gl.getResultsArray(); // java map array
		//if (i == (slices-1)) gl.writetoResultsTable(true);
		//else gl.writetoResultsTable(false);
		IJ.showProgress((i+1)/slices);
		status = "" + i + "/"+  slices;
		IJ.showStatus(status);
	}
	IJ.showStatus("GLCM calculated: writing results");

	var rt = ResultsTable.getResultsTable();
	rt.reset();
	for (var i = 0; i < slices; i++){
		/*
		gl.setglcm(glcmA[i]);
		if (i == (slices-1)) gl.writetoResultsTable(true);
		else gl.writetoResultsTable(false);
		*/
		WritetoResults(pA[i], rt, paranamesA)
		status = "" + i + "/"+  slices;
		IJ.showStatus("GLCM calculated: writing results:" + status);		
		IJ.showProgress((i+1)/slices);
	}
	rt.show("Results");
	IJ.showStatus("GLCM Done");
	
}

//map is an instance of javamap, retrieved from the java instance
// rt is the current results table. 
//paranamesA is the arary containing parameter names
function WritetoResults(map, rt, paranamesA){
	var javamap = HashMap(map);
	//IJ.log(javamap.get("Energy"));
	var row = rt.getCounter();
	rt.incrementCounter();
	for (var k = 0; k < paranamesA.length; k++)
		rt.setValue(paranamesA[k], row, javamap.get(paranamesA[k]));	
}
