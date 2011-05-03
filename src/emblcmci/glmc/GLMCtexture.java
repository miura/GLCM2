package emblcmci.glmc;
/** Kota Miura miura@embl.de
 * 		most of the parts copied from GLCM_TextToo.java written by Toby C. Cornish @jhmi.edu
 * 		see GLCM_TextToo.java header part for details. 
 * 	for use from scripts and other class as stand alone:
 * 		use constructor GLMCtexture(int d, int phi, boolean  symmetry)
 * 		use calcGLMC(ImagePlus imp) to set imageplus object and calculate GLMC
 */

/* sample javascript 

importClass(Packages.emblcmci.glmc.GLMCtexture);
g = new GLMCtexture(IJ.getImage(), 2, 45, true, true);
glmc = g.calcGLMC();
ht = g.getResultsArray(g);
IJ.log(ht.get("Contrast")) ;

 */

import java.awt.Rectangle;
import java.util.HashMap;
import java.util.Map;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.process.ImageProcessor;

public class GLMCtexture {
	static int d = 1;
	static int phi = 0;
	static boolean symmetry = true;
	static boolean doASM = true;
	static boolean doContrast = true;
	static boolean doCorrelation = true;
	static boolean doIDM = true;
	static boolean doEntropy = true;
	static boolean doEnergy = true;
	static boolean doInertia = true;
	static boolean doHomogeneity = true;
	static boolean doProminence = true;
	static boolean doVariance = true;
	static boolean doShade = true;
	
	boolean rt_reset = true;


	ResultsTable rt = ResultsTable.getResultsTable();

	double [][] glcm; 
	double meanx=0.0;
	double meany=0.0;
	double stdevx=0.0;
	double stdevy=0.0;

	public double Asm;
	public double Idm;
	public double Contrast;
	public double Energy;
	public double Entropy;
	public double Homogeneity;
	public double Variance;
	public double Shade;
	public double Prominence;
	public double Inertia;
	public double Correlation;
	public double GLCMsum;
	
	ImagePlus imp;

	/**Constructor for use as library
	 * all the parameters will be by default be measured (true). 
	 * @param d
	 * @param phi
	 * @param symmetry
	 */
	@SuppressWarnings("static-access")
	public GLMCtexture(int d, int phi, boolean  symmetry, boolean rt_reset){
		this.d = d;
		this.phi = phi;
		this.symmetry = symmetry;
		this.rt_reset = rt_reset;
	}
	
	@SuppressWarnings("static-access")
	public GLMCtexture(ImagePlus imp, int d, int phi, boolean  symmetry, boolean rt_reset){
		this.imp = imp;
		this.d = d;
		this.phi = phi;
		this.symmetry = symmetry;
		this.rt_reset = rt_reset;
	}
	
	public GLMCtexture() {
		// TODO Auto-generated constructor stub
	}

	public static void setD(int d) {
		GLMCtexture.d = d;
	}

	public static void setPhi(int phi) {
		GLMCtexture.phi = phi;
	}	

	public void setRt_reset(boolean rtReset) {
		rt_reset = rtReset;
	}
	
	public void doBasicStats(){
		double [] px = new double [256];
		double [] py = new double [256];
		meanx=0.0;
		meany=0.0;
		stdevx=0.0;
		stdevy=0.0;

		// Px(i) and Py(j) are the marginal-probability matrix; sum rows (px) or columns (py) 
		// First, initialize the arrays to 0
		for (int i=0;  i<256; i++){
			px[i] = 0.0;
			py[i] = 0.0;
		}

		// sum the glcm rows to Px(i)
		for (int i=0;  i<256; i++) {
			for (int j=0; j<256; j++) {
				px[i] += glcm [i][j];
			} 
		}

		// sum the glcm rows to Py(j)
		for (int j=0;  j<256; j++) {
			for (int i=0; i<256; i++) {
				py[j] += glcm [i][j];
			} 
		}

		// calculate meanx and meany
		for (int i=0;  i<256; i++) {
			meanx += (i*px[i]);
			meany += (i*py[i]);
		}

		// calculate stdevx and stdevy
		for (int i=0;  i<256; i++) {
			stdevx += ((Math.pow((i-meanx),2))*px[i]);
			stdevy += ((Math.pow((i-meany),2))*py[i]);
		}
	}
		
	//=====================================================================================================
	// This is the generic moments function from parker -- may implement in the future
	// k is the power for the moment

	/*
if (doMoments == true){
	double y=0.0;
	double z;
	double k;
	double moments = 0.0;

	for (int i=0;  i<256; i++)  {
		for (int j=0; j<256; j++) {
			if (k>0) {
				z = Math.pow ((i-j), k);
			} else {
				if (i == j) continue;
				z = Math.pow ((i-j), -1*k);
				z = glcm[i][j]/z;
			}
			moments += z * glcm[i][j];
		}
	}
	rt.setValue("Angular Second Moment", row, asm);
}
	 */	
	//===============================================================================================
	// calculate the inverse difference moment (idm) (Walker, et al. 1995)
	// this is calculated using the same formula as Conners, et al., 1984 "Local Homogeneity"

	//=====================================================================================================
	// calculate the angular second moment (asm)
	
	public double getAngular2ndMoment(){
		double asm = 0.0;
		for (int i=0;  i<256; i++)  {
			for (int j=0; j<256; j++) {
				asm += (glcm[i][j]*glcm[i][j]);
			}
		}
		return asm;
	}

	public double getIDM(){
		double IDM = 0.0;
		for (int i=0;  i<256; i++)  {
			for (int j=0; j<256; j++) {
				IDM += ((1/(1+(Math.pow(i-j,2))))*glcm[i][j]);
			}
		}
		return IDM;
	}
	//=====================================================================================================
	// calculate the angular second moment (asm)
	
	public double getContrast(){
		double contrast=0.0;

		for (int i=0;  i<256; i++)  {
			for (int j=0; j<256; j++) {
				contrast += Math.pow(Math.abs(i-j),2)*(glcm[i][j]);
			}
		}
		return contrast;
	}
	//===============================================================================================
	// calculate the energy	
	public double getEnergy(){
		double energy = 0.0;
		for (int i=0;  i<256; i++)  {
			for (int j=0; j<256; j++) {
				energy += Math.pow(glcm[i][j],2);
			}
		}
		return energy;
	}
	//===============================================================================================
	// calculate the entropy (Haralick et al., 1973; Walker, et al., 1995)
	public double getEntropy(){
		double entropy = 0.0;
		for (int i=0;  i<256; i++)  {
			for (int j=0; j<256; j++) {
				if (glcm[i][j] != 0) {
					entropy = entropy-(glcm[i][j]*(Math.log(glcm[i][j])));
					//the next line is how Xite calculates it -- I am not sure why they use this, I do not think it is correct
					//(they also use log base 10, which I need to implement)
					//entropy = entropy-(glcm[i][j]*((Math.log(glcm[i][j]))/Math.log(2.0)) );
				}
			}
		}
		return entropy;
	}
	//===============================================================================================
	// calculate the homogeneity (Parker)
	// "Local Homogeneity" from Conners, et al., 1984 is calculated the same as IDM above
	// Parker's implementation is below; absolute value of i-j is taken rather than square
	
	public double getHomogeneity(){
		double homogeneity = 0.0;
		for (int i=0;  i<256; i++) {
			for (int j=0; j<256; j++) {
				homogeneity += glcm[i][j]/(1.0+Math.abs(i-j));
			}
		}
		return homogeneity;
	}

	//===============================================================================================
	// calculate the variance ("variance" in Walker 1995; "Sum of Squares: Variance" in Haralick 1973)

	public double getVariance(){
		double variance = 0.0;
		double mean = 0.0;

		mean = (meanx + meany)/2;
		/*
	// this is based on xite, and is much greater than the actual mean -- it is here for reference only
	for (int i=0;  i<256; i++)  {
		for (int j=0; j<256; j++) {
			mean += glcm[i][j]*i*j;
		}
	}
		 */

		for (int i=0;  i<256; i++)  {
			for (int j=0; j<256; j++) {
				variance += (Math.pow((i-mean),2)* glcm[i][j]);
			}
		}
		return variance;		
	}
	
	//===============================================================================================
	// calculate the shade (Walker, et al., 1995; Connors, et al. 1984)
	public double getShade(){
		double shade = 0.0;

		// calculate the shade parameter
		for (int i=0;  i<256; i++) {
			for (int j=0; j<256; j++) {
				shade += (Math.pow((i+j-meanx-meany),3)*glcm[i][j]);
			}
		}
		return shade;
	}

	//==============================================================================================
	// calculate the prominence (Walker, et al., 1995; Connors, et al. 1984)	
	public double getProminence(){
		double prominence=0.0;

		for (int i=0;  i<256; i++) {
			for (int j=0; j<256; j++) {
				prominence += (Math.pow((i+j-meanx-meany),4)*glcm[i][j]);
			}
		}
		return prominence;
	}
	
	//===============================================================================================
	// calculate the inertia (Walker, et al., 1995; Connors, et al. 1984)	
	public double getInertia(){
		double inertia = 0.0;
		for (int i=0;  i<256; i++)  {
			for (int j=0; j<256; j++) {
				if (glcm[i][j] != 0) {
					inertia += (Math.pow((i-j),2)*glcm[i][j]);
				}
			}
		}
		return inertia;
	}
	//=====================================================================================================
	// calculate the correlation
	// methods based on Haralick 1973 (and MatLab), Walker 1995 are included below
	// Haralick/Matlab result reported for correlation currently; will give Walker as an option in the future

	public double getCorrelation(){
		double correlation=0.0;

		// calculate the correlation parameter
		for (int i=0;  i<256; i++) {
			for (int j=0; j<256; j++) {
				//Walker, et al. 1995 (matches Xite)
				//correlation += ((((i-meanx)*(j-meany))/Math.sqrt(stdevx*stdevy))*glcm[i][j]);
				//Haralick, et al. 1973 (continued below outside loop; matches original GLCM_Texture)
				//correlation += (i*j)*glcm[i][j];
				//matlab's rephrasing of Haralick 1973; produces the same result as Haralick 1973
				correlation += ((((i-meanx)*(j-meany))/( stdevx*stdevy))*glcm[i][j]);
			}
		}
		//Haralick, et al. 1973, original method continued.
		//correlation = (correlation -(meanx*meany))/(stdevx*stdevy);
		return correlation;
	}

	//===============================================================================================
	// calculate the sum of all glcm elements
	public double getGLCMsum(){
		double sum = 0.0;
		for (int i=0; i<256; i++)  {
			for (int j=0; j<256; j++) {
				sum = sum + glcm[i][j];
			}
		}
		return sum;
	}
	//=========================================================================================

	public void setFieldParameters(){
		this.Asm = getAngular2ndMoment();
		this.Idm = getIDM();
		this.Contrast = getContrast();
		this.Energy = getEnergy();
		this.Entropy = getEntropy();
		this.Homogeneity = getHomogeneity();
		this.Variance = getVariance();
		this.Shade = getShade();
		this.Prominence = getProminence();
		this.Inertia = getInertia();
		this.Correlation = getCorrelation();
		this.GLCMsum = getGLCMsum();		
	}
	
	// implementation of the dialog
	public boolean showDialog() {
		GenericDialog gd = new GenericDialog("GLCM Texture v0.001");
		gd.addNumericField ("Enter the size of the step in pixels",  d, 0);

		String [] angles={"0", "45", "90", "135"};
		gd.addChoice("Select the direction of the step", angles, Integer.toString(phi));
		gd.addCheckbox("Symmetrical GLCM?", symmetry);
		gd.addMessage("Calculate which parameters?");   
		gd.addCheckbox("Angular Second Moment  ", doASM);
		gd.addCheckbox("Contrast  ", doContrast);
		gd.addCheckbox ("Correlation  ", doCorrelation);
		gd.addCheckbox ("Inverse Difference Moment  ", doIDM);
		gd.addCheckbox ("Entropy   ", doEntropy);
		gd.addCheckbox ("Energy   ", doEnergy);
		gd.addCheckbox ("Inertia   ", doInertia);
		gd.addCheckbox ("Homogeneity   ", doHomogeneity);
		gd.addCheckbox ("Prominence   ", doProminence);
		gd.addCheckbox ("Variance   ", doVariance);
		gd.addCheckbox ("Shade   ", doShade);

		gd.showDialog();
		if (gd.wasCanceled()) return false;

		d=(int) gd.getNextNumber();
		phi=Integer.parseInt(gd.getNextChoice());
		symmetry=gd.getNextBoolean();
		doASM=gd.getNextBoolean();
		doContrast=gd.getNextBoolean();
		doCorrelation=gd.getNextBoolean();
		doIDM=gd.getNextBoolean();
		doEntropy=gd.getNextBoolean();
		doEnergy=gd.getNextBoolean();
		doInertia=gd.getNextBoolean();
		doHomogeneity=gd.getNextBoolean();
		doProminence=gd.getNextBoolean();
		doVariance=gd.getNextBoolean();
		doShade=gd.getNextBoolean();

		return true;
	}

	public double [][] calcGLMC(){
		if (imp == null)
			imp = IJ.getImage();
		return calcGLMC(imp);
	}
	
	public double [][] calcGLMC(ImagePlus imp){
		return calcGLMC(imp.getProcessor());
	}
	
	/**main part that does the calculation of GLMC
	 * 
	 * @param ip
	 */
	public double [][] calcGLMC(ImageProcessor ip){
		// use the bounding rectangle ROI to roughly limit processing
		Rectangle roi = ip.getRoi();

		// get byte arrays for the image pixels and mask pixels
		int width = ip.getWidth();
		int height = ip.getHeight();
		byte [] pixels = (byte []) ip.getPixels();
		byte [] mask = ip.getMaskArray();

		// value = value at pixel of interest; dValue = value of pixel at offset    
		int value;
		int dValue;
		double totalPixels = roi.height * roi.width;
		if (symmetry) totalPixels = totalPixels * 2; 
		double pixelProgress = 0;
		double pixelCount = 0;

		//====================================================================================================
		// compute the Gray Level Correlation Matrix

		int offsetX = 1;
		int offsetY = 0;
		glcm = new double [256][256];

		// set our offsets based on the selected angle
//		if (phi == 0) {
//			offsetX = d;
//			offsetY = 0;
//		} else if (phi == 45) {
//			offsetX = d;
//			offsetY = -d;
//		} else if (phi == 90) {
//			offsetX = 0;
//			offsetY = -d;
//		} else if (phi == 135) {
//			offsetX = -d;
//			offsetY = -d;
//		} else {
//			// the angle is not one of the options
//			IJ.showMessage("The requested angle,"+phi+", is not one of the supported angles (0,45,90,135)");
//		}
		double rad = Math.toRadians(-1.0 * phi);
		offsetX = (int) ((int) d* Math.round(Math.cos(rad)));
		offsetY = (int) ((int) d* Math.round(Math.sin(rad)));
		

		// loop through the pixels in the ROI bounding rectangle
		for (int y=roi.y; y<(roi.y + roi.height); y++) 	{
			for (int x=roi.x; x<(roi.x + roi.width); x++)	 {
				// check to see if the pixel is in the mask (if it exists)
				if ((mask == null) || ((0xff & mask[(((y-roi.y)*roi.width)+(x-roi.x))]) > 0) ) {
					// check to see if the offset pixel is in the roi
					int dx = x + offsetX;
					int dy = y + offsetY;
					if ( ((dx >= roi.x) && (dx < (roi.x+roi.width))) && ((dy >= roi.y) && (dy < (roi.y+roi.height))) ) {
						// check to see if the offset pixel is in the mask (if it exists) 
						if ((mask == null) || ((0xff & mask[(((dy-roi.y)*roi.width)+(dx-roi.x))]) > 0) ) {
							value = 0xff & pixels[(y*width)+x];
							dValue = 0xff & pixels[(dy*width) + dx];
							glcm [value][dValue]++;		  			
							pixelCount++;
						}
						// if symmetry is selected, invert the offsets and go through the process again
						if (symmetry) {
							dx = x - offsetX;
							dy = y - offsetY;
							if ( ((dx >= roi.x) && (dx < (roi.x+roi.width))) && ((dy >= roi.y) && (dy < (roi.y+roi.height))) ) {
								// check to see if the offset pixel is in the mask (if it exists) 
								if ((mask == null) || ((0xff & mask[(((dy-roi.y)*roi.width)+(dx-roi.x))]) > 0) ) {
									value = 0xff & pixels[(y*width)+x];
									dValue = 0xff & pixels[(dy*width) + dx];
									glcm [dValue][value]++;		  			
									pixelCount++;
								}	
							}
						}
					}  
				}
				pixelProgress++;	
				IJ.showProgress(pixelProgress/totalPixels);
			}
		}

		// convert the GLCM from absolute counts to probabilities
		for (int i=0; i<256; i++)  {
			for (int j=0; j<256; j++) {
				glcm[i][j] = (glcm[i][j])/(pixelCount);
			}
		}
		return glcm;
	}
	
	public void writetoResultsTable(GLMCtexture gl){
		ResultsTable rt = ResultsTable.getResultsTable();
		if (this.rt_reset) rt.reset();
		gl.doBasicStats();
		gl.setFieldParameters();
		int row = rt.getCounter();	
		rt.incrementCounter();
		if (doASM)			
			rt.setValue("Angular Second Moment", row, gl.getAngular2ndMoment());
		if (doIDM)
			rt.setValue("Inverse Difference Moment", row, gl.getIDM());
		if (doContrast)
			rt.setValue("Contrast", row, gl.getContrast());	
		if (doEnergy)
			rt.setValue("Energy", row, gl.getEnergy());
		if (doEntropy)
			rt.setValue("Entropy", row, gl.getEntropy());
		if (doHomogeneity)
			rt.setValue("Homogeneity", row, gl.getHomogeneity());		
		if (doVariance)
			rt.setValue("Variance", row, gl.getVariance());
		if (doShade)
			rt.setValue("Shade", row, gl.getShade());
		if (doProminence)
			rt.setValue("Prominence", row, gl.getProminence()); 
		if (doInertia)
			rt.setValue("Inertia", row, gl.getInertia());
		if (doCorrelation)
			rt.setValue("Correlation", row, gl.getCorrelation()); 
		
		rt.setValue("Sum of all GLCM elements", row, gl.getGLCMsum());
		rt.show("Results");		
	}
	
	public HashMap<?, ?> getResultsArray(GLMCtexture gl){
		HashMap<String, Double> res = new HashMap<String, Double>();
		res.put("Angular Second Moment", gl.getAngular2ndMoment());
		res.put("Inverse Difference Moment", gl.getIDM());
		res.put("Contrast", gl.getContrast());
		res.put("Energy", gl.getEnergy());
		res.put("Homogeneity", gl.getHomogeneity());
		res.put("Variance", gl.getVariance());
		res.put("Shade", gl.getShade());
		res.put("Prominence", gl.getProminence());
		res.put("Inertia", gl.getInertia());
		res.put("Correlation", gl.getCorrelation());
		res.put("Sum of all GLCM elements", gl.getGLCMsum());
		return res;
	}
	

}
