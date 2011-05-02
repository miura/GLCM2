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

/* 5/2/2011 Modified by Kota Miura (miura@embl.de) 
 * 	added setters, getters
 * 	out-sourced  processes in run() to several methods. 
 * 	major part moved to a new class emblcmci.GLMC.GLMCtexture to enable library-like usage
 */

import ij.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import emblcmci.glmc.GLMCtexture;

//==========================================================
public class GLCM_TextureToo implements PlugInFilter {
	
	GLMCtexture gl = new GLMCtexture();
	
	public int setup(String arg, ImagePlus imp) {
		if (imp!=null && !gl.showDialog()) return DONE;
		//perhaps not reseting the resultsTable would be better... ??
		return DOES_8G+DOES_STACKS+SUPPORTS_MASKING;
	}

	public void run(ImageProcessor ip) {
		//calculates GLMC
		gl.calcGLMC(ip);
		// show results in the Results Table
		gl.writetoResultsTable(gl, true);
	}
}
