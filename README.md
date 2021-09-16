# GLCM2

Calculates GLCM texture parameters. 
The original plugin [GLCM Texture](https://imagej.nih.gov/ij/plugins/texture.html) was originally written by Julio E. Cabrera (jcabrera at mail.nih.gov). 

With this plugin, the original plugin was 
1. Mavenized. 
2. Allows a better access from scripting. 

## Other similar ImageJ plugins (as of Sept. 2021)

1. [pyRadiomics_GUI](https://github.com/salimkanoun/pyRadiomics_GUI)
  - this is a wrapper for a python library, to use it as an ImageJ plugin. 
  - If you prefere to use it directly in python, see here https://pyradiomics.readthedocs.io/en/latest/installation.html
2. [ImageJ-OPs has implementations](https://javadoc.scijava.org/ImageJ/net/imagej/ops/features/haralick/package-summary.html). 
  - See [the forum thread](https://forum.image.sc/t/glcm-texturetool-latest-version/9386)
  - There is a very good discussion on how to use these Ops implementations. It's about using it for CLIJ, but it actually is a very general discussion. (CLIJ optimisation on GLCM - Image Analysis](https://forum.image.sc/t/clij-optimisation-on-glcm/41240/2]
