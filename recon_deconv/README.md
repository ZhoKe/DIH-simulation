# Deconvolution with Inverse Iterative Particle Extraction matlab
Matlab version of the deconvolution method developed by Mostafa Toloui for DIH-PTV in densely-seeded flows. May be applied to other applications as well providing the particles are point-like, uniformly spatially distributed, and imaged at a relatively dense concentration. The precise bounds are not currently known. The size of the image that can be processed is limited by computer RAM, with a maximum size of approximately 512x512x512.

# Processing Steps:
1.	SKSMT_IIPE_v08.m

# Inputs:
All inputs parameters are defined at the beginning of SKSMT_IIPE_v08.m. they are summarized here to define them and explain their use.
## Primary Parameters
- org: Flag to determine whether to perform a step using the standard (original) reconstruction method. Default is 0 (false). If set to true, will output reconstructed planes and projected images and terminate without proceeding further. Useful for checking reconstruction parameters.
- perform_background_subtraction: Flag to identify whether input images need to be enhanced. If set to true, the inputs in Raw_img will be enhanced by subtracting their mean background and saved to Enh_img. If set to false, the images described by Enh_img will be used directly. Note that the background subtraction will only work well if the Sequence of images specified is sufficiently long to compute a reasonable mean background.
- Niteration: Number of IIPE iterations to be performed. Default is 6. Dense particle fields may require more while sparse fields may require fewer.
- RI: Refractive index of the fluid (i.e. 1.33 for water)
- Dp: Estimated particle diameter (units: pixels)
- Lambda_air: Wavelength of imaging laser light in air (i.e. the normal wavelength of the lasser)
- Reso: Image resolution (i.e. pixel size) (units: nm/pixel)
- ddz: Reconstruction plane spacing (units: nm)
- Z1: Depth offset of first reconstructed plane (units: nm)
- nz: Number of reconstruction planes to use
- Sequence: List of image frame numbers to be processed
- ImgSiz: Size of the input image in pixels. [rows, columns]
- Beta: Deconvolution relaxation parameter to avoid division by zero. Default of 0.1 rarely needs to be adjusted
- pathn: File path to the project directory. Output folders will be created here.
- Raw_img: String defining the format of the raw hologram image file name. Has one format component to specify the image number in the sequence.
- Enh_img: String defining the format of the enhanced hologram image file name. Must contain two format components (e.g. %05d). The first specifies the image number in the sequence (corresponding to the Sequence parameter). The second is the iteration number, the input image stats with 1. It is generally assumed that the input images will be located in the pathn/Holograms folder.
## Infrequently Modified Parameters
- R_dial: Radius of dilation disk used in particle removal. Default of 2 rarely adjusted
- C0: Scaling factor to apply to mask replacement intensity value. Default of 1 rarely adjusted
- Niter: Number of times to remove each particle from the hologram. Default of 4 rarely adjusted
- Corast_recovery: Flag to indicate whether contrast recovery should be performed in the removal step. Nearly always set to true.
- FA: Flag indicating modification to removal using additional dilation. Nearly always set to false
- C0_STA: Enhancement parameter. Scale to the mean of the standard deviation, influences spread of the output histogram.
- LL_GS: Enhancement parameter. Lowest level grayscale shift applied to output image. 

# Outputs:
The following output folders are created and filled during code execution: 3D_imgs, Cmb_imgs, Dv_rec, and (optional) Org_rec. Additionally, some outputs are placed in the Holograms folder alongside the input images.

The primary output of the code is the particle centroid list located in **3D_imgs/ParticleCentroid_%03d.txt**. This is a space-delimited text file containing the x, y, and z centroid coordinates for each particle along with the bounding box. There are two header lines identifying the data columns.

All other outputs are intermediate steps. These files are needed in the processing. However, due to the large number of small files generated they may use a substantial amount of storage space. It is suggested that these outputs be deleted after processing if they are not going to be used further.
