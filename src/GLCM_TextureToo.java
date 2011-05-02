//=================================================================================================
// GLCM_Texture_Too, v. 0.008
// Toby C. Cornish, tcornis3@jhmi.edu (or toby@tobycornish.com)
// 11/26/2007
// modified from GLCM_Texture (Gray Level Correlation Matrix Texture Analyzer) v0.4, Author: Julio E. Cabrera, 06/10/05
// 
// CHANGELOG:
// ------------------------------------------------------------------------------------------------
// 11/26/07 GLCM_TextureToo v0.008 - minor fixes; added references section -- tc
// 11/22/07 GLCM_TextureToo v0.007 - minor fixes to parameters -- tc
//                                 - moved mean and stdev calculations to common area -- tc
// 11/20/07 GLCM_TextureToo v0.006 - confirmed results against parker and Xite -- tc
//                                 - added preliminary shade, prominence, variance, homogeneity, inertia
//                                 - differentiated Haralick correlation from Walker correlation
// 11/19/07 GLCM_TextureToo v0.005 - changed method of calculating correlation -- tc
//                                 - should be closest in nomenclature to Walker, et al.  1995
// 11/19/07 GLCM_TextureToo v0.004 - corrected column name of idm -- tc
// 11/13/07 GLCM_TextureToo v0.003 - changed from roi.contains() to byte[] checking -> much faster -- tc
// 11/13/07 GLCM_TextureToo v0.002 - added progress bar --tc
// 11/11/07 GLCM_TextureToo v0.001 - inherited portions of the codebase from GLCM_Texture 0.4 -- tc
//                                 - fundamental rewrite of GLCM calculation
//                                 - added irregular ROI support
//                                 - added symmetrical/non-symmetrical GLCM calculations
//                                 - corrected directionality (phi) to be 0,45,90,135
//
//=================================================================================================
//
// References: 
//   R.M. Haralick, Texture feature for image classification, IEEE Trans. SMC 3 (1973) (1), pp. 610–621.
//   Conners, R.W., Trivedi, M.M., and Harlow, C.A., Segmentation of a High-Resolution Urban Scene
//     Using Texture Operators, CVGIP(25), No. 3, March, 1984, pp. 273-310.
//   Walker, RF, Jackway, P and Longstaff, ID (1995) Improving Co-occurrence Matrix Feature Discrimination.'
//     In  DICTA '95, 3rd Conference on Digital Image Computing: Techniques and Application, 6 - 8 December,
//     1995, pages 643-648.
//   Parker, JR, Algorithms for Image Processing and Computer Vision, John Wiley & Sons, 1997.
//   Image processing lab, Department of Informatics, University of Oslo. Xite v1.35: glcmParameter.c, v1.30
//     2004/05/05 07:34:19 (2004)

/* Modified by Kota Miura (miura@embl.de)
 * setters, getters
 * created several method to outsource. 
 * 
 */

import ij.*;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.awt.*;

import emblcmci.glmc.GLMCtexture;
import ij.measure.ResultsTable;

//==========================================================
public class GLCM_TextureToo implements PlugInFilter {
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
	
	GLMCtexture glmc = new GLMCtexture();
	
	public int setup(String arg, ImagePlus imp) {
		if (imp!=null && !glmc.showDialog()) return DONE;
		//perhaps not reseting the resultsTable would be better... ??
		rt.reset();
		return DOES_8G+DOES_STACKS+SUPPORTS_MASKING;
	}

	public void run(ImageProcessor ip) {

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


		//=====================================================================================================

		// convert the GLCM from absolute counts to probabilities
		for (int i=0; i<256; i++)  {
			for (int j=0; j<256; j++) {
				glcm[i][j] = (glcm[i][j])/(pixelCount);
			}
		}
		glmc.doBasicStats();
		glmc.setFieldParameters();
		int row = rt.getCounter();	
		rt.incrementCounter();
		if (doASM)			
			rt.setValue("Angular Second Moment", row, glmc.getAngular2ndMoment());
		if (doIDM)
			rt.setValue("Inverse Difference Moment", row, glmc.getIDM());
		if (doContrast)
			rt.setValue("Contrast", row, glmc.getContrast());	
		if (doEnergy)
			rt.setValue("Energy", row, glmc.getEnergy());
		if (doEntropy)
			rt.setValue("Entropy", row, glmc.getEntropy());
		if (doHomogeneity)
			rt.setValue("Homogeneity", row, glmc.getHomogeneity());		
		if (doVariance)
			rt.setValue("Variance", row, glmc.getVariance());
		if (doShade)
			rt.setValue("Shade", row, glmc.getShade());
		if (doProminence)
			rt.setValue("Prominence", row, glmc.getProminence()); 
		if (doInertia)
			rt.setValue("Inertia", row, glmc.getInertia());
		if (doCorrelation)
			rt.setValue("Correlation", row, glmc.getCorrelation()); 
		
		rt.setValue("Sum of all GLCM elements", row, glmc.getGLCMsum());
		rt.show("Results");
	}


}
