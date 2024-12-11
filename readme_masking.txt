# Creating Mask of AOI for PSI processing
#
# Author: Mahmut Arikan <m.arikan@tudelft.nl>
# Date:	  20060711
# Version: 0.11
# Updated: 20071106, FvL

#General outline of creating mask
1. Open multireflectivity map in a Raster Editor (Gimp or Grass)
	- Important: Don't use flipped multireflectivity map images!
2. Add a layer and create a mask using brush tools on the area of
interest (AOI). Use white brush to draw AOI and then  flood fill the
rest with black color.  
3. Save resulting mask to binary mask and then convert it to
PS_tools_box compatible .raw file (format unsigned integer*8)

# Example using gimp
0. Create a mrm.ras file from the mrm.raw file using cpxfiddle:
   cpxfiddle -w[width] -fr4 -M1/1 -qnormal -osunraster -cgray mrm.raw > mrm.ras
1. Open mrm.ras in GIMP 
2. Enhance contrast (optional step)
	- MENU: LAYER/COLORS/LEVELS , just include main part of histogram
3. Add new transparent layer 
	- MENU: LAYER/NEW LAYER
	- if layer control window is not visible activate it using [CRTL+L]
	- lock background layer by click area in between eye symbol and layer
	  icon in layer control window (chain icon will appear when you click).
4. Using pencil brush draw over and around area interest with white color.
	- Create a new brush if size of brush is small
	  - In Brushes window on right-click menu select new brush and create a
	    new one with size 76.
	- If necessary, change opacity of the layer after you have
	created mask.
	- If necessary, you can flood fill the area of interest with
	white color.

5. Remove background layer
	- Select background from layer control window and click trash icon.
6. After you have completed drawing flood fill rest of the mask with black color
	- MENU: SELECT/BY COLOR and click in transparant region
	- MENU: EDIT/FILL WITH BG COLOR
7. Flatten created mask layer
	- MENU: IMAGE/FLATTEN iMAGE
8. Switched to color mode indexed
	- MENU: IMAGE/MODE/INDEXED and select Use black and white (1-bit) palette
9. Save as tif or any other format compatible with matlab

# Converting mask image to .raw file (unsigned integer*8)
1. In matlab use im2PSraw <input.tif> <output.raw>
	- im2PSraw mask.tif mask.raw 
2. You are done and just copy mask.raw under your PSI processing directory
	- You can check mask.raw using any raw image reader like OpenEV 

#EOF
