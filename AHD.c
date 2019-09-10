/*
 * AHD.c
 * 
 * Copyright 2019 qingyao <qingyao@blink-7572-1>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#define ushort UshORt
typedef unsigned char uchar;
typedef unsigned short ushort;

#define TS 256		/* Tile Size */
#define height 5472
#define width 7296 
#define colors 3
#define filters 0x94949494  //color_desc =  b'RGBG'
ushort (*image)[colors];
ushort *raw_img;

ushort (*pix)[colors], (*rix)[3]; // pix => original raw data ; rix => after interpolation 
int i, j, k, top, left, row, col, tr, tc, c, d, z, val, hm[2],color; 
static const int dir[4] = { -1, 1, -TS, TS }; // direction => {left, right, up, down}
unsigned l_diff[2][4], ab_diff[2][4], l_eps, ab_eps; // L difference, color difference and epsilon (parameters)
float r, xyz[3], xyz_cam[3][colors], rgb_cam[3][colors]; 
ushort (*rgb)[TS][TS][3];
short (*lab)[TS][TS][3], (*lix)[3];
char (*homo)[TS][TS], *buffer;
const double xyz_rgb[3][3] = {			 // XYZ from RGB 
    { 0.412453, 0.357580, 0.180423 },    // RGB Working Space : sRGB 
    { 0.212671, 0.715160, 0.072169 },    // Reference White : D65 
    { 0.019334, 0.119193, 0.950227 } };  
double d65[3] = { 0.950456, 1, 1.088754 }; 
/*int filter[8] = {0x16161616,0x61616161,0x49494949,0x94949494,
	               0x1e1e1e1e,0xe1e1e1e1,0x4b4b4b4b,0xb4b4b4b4};
  All RGB cameras use one of these Bayer grids :
    if filters = 0x16161616:	  0x61616161:	  0x49494949:	 0x94949494:  (3 colors)
                 0 1 2 3 4 5	  0 1 2 3 4 5	  0 1 2 3 4 5	  0 1 2 3 4 5
	           0 B G B G B G	0 G R G R G R	0 G B G B G B	0 R G R G R G
	           1 G R G R G R	1 B G B G B G	1 R G R G R G	1 G B G B G B
	           2 B G B G B G	2 G R G R G R	2 G B G B G B	2 R G R G R G
	           3 G R G R G R	3 B G B G B G	3 R G R G R G	3 G B G B G B
	            
	if filters = 0x1e1e1e1e:	  0xe1e1e1e1:	  0x4b4b4b4b:	 0xb4b4b4b4:  (4 colors)
                 0 1 2 3 4 5	  0 1 2 3 4 5	  0 1 2 3 4 5	  0 1 2 3 4 5 
	           0 B g B g B g	0 G R G R G R	0 g B g B g B	0 R G R G R G
	           1 G R G R G R	1 B g B g B g	1 R G R G R G	1 g B g B g B
	           2 B g B g B g	2 G R G R G R	2 g B g B g B	2 R G R G R G
	           3 G R G R G R	3 B g B g B g	3 R G R G R G	3 g B g B g B
*/
// Define Functions
#define FORC3 for (c=0; c < 3; c++)
#define FORC4 for (c=0; c < 4; c++)
#define FORCC for (c=0; c < colors; c++)
#define SQR(x) ((x)*(x))
#define ABS(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define LIM(x,min,max) MAX(min,MIN(x,max)) 
#define CLIP(x) LIM(x,0,4095)   // up bound = 4095 for 12 bit image
#define FC(row,col)  (filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)  // return color index of selected pixel   
#define BAYER(row,col)  image[row * width + col][FC(row,col)] // return known raw data of selected pixel

// Border interpolation before applying AHD
void border_interpolate (unsigned int border)
{
	unsigned row, col, y, x, c, sum[8];
	for (row=0; row < height; row++){
		for (col=0; col < width; col++){
			if (col==border && row >= border && row < height-border){
				col = width-border;
			}
		    memset (sum, 0, sizeof(sum)); // initialize sum with all 0
		    for (x=row-1; x != row+2; x++){
				for (y=col-1; y != col+2; y++){
					if (x < height && y < width){
						sum[FC(x,y)] += BAYER(x,y); 
						sum[FC(x,y)+4]++;
					}}}
			FORCC {if (c != FC(row,col)){
				image[row*width+col][c] = sum[c] / sum[c+4];
			}}}}
}

// ADAPTIVE HOMOGENEITY-DIRECTED (AHD) DEMOSAICING ALGORITHM 
void ahd_demosaic()
{
	float cbrt[0x10000];
	
	for (i=0; i < 0x10000; i++) {
		r = i / 65535.0;
		cbrt[i] = r > 0.008856 ? pow((double)r,1/3.0) : 7.787*r + 16/116.0;
	}
	
	static const float table[12] = 

	{ -1.936280,  1.800443, -1.448486,  2.584324,
	   1.405365, -0.524955, -0.289090,  0.408680,
	  -1.204965,  1.082304,  2.941367, -1.818705 };

	for (i=0; i < 3; i++){
		FORCC {rgb_cam[i][c] = table[i*4+c];}}
	
	// rgb_cam : Color trans matrix to transform pixel value from camera color space to srgb color space.

	for (i=0; i < 3; i++){  // cam to xyz coefficient
		for (j=0; j < colors; j++){
			for (xyz_cam[i][j] = k=0; k < 3; k++){
				xyz_cam[i][j] += xyz_rgb[i][k] * rgb_cam[k][j] / d65[i];
			}
		}
	}
	border_interpolate (6);
	buffer = (char *) malloc (26*TS*TS);		/* Assign 1664 kB and return the address of first byte*/
	rgb  = (ushort(*)[TS][TS][3]) buffer;            // RGB img
	lab  = (short(*)[TS][TS][3])(buffer + 12*TS*TS);// LAB img = buffer + 2 x 2 x 3 x TS x TS, stored after RGB
	                                                // 2 directions, 3 color channels, size TS x TS, each point 2 Bytes
	homo = (char  (*)[TS][TS])   (buffer + 24*TS*TS);// homo map
	printf("Start to Demosaic.\n");
	for (top=3; top < height-6; top += TS-7){        // Each processing image pattern has size 249x249			   
		for (left=3; left < width-6; left += TS-7){  // The first pattern starts with (3,3) to (249,249)						
			memset (rgb, 0, 12*TS*TS);
			/*  Interpolate green horizontally and vertically:		*/
			for (row = top; row < top+TS && row < height-3; row++){
				col = left + (FC(row,left) & 1); // locate B and R for interpolating G on these positions 
				/* (3,3)            FC(row,col) = 2              B   
				           (4,4)    FC(row,col) =   0                R
                   (5,3)            FC(row,col) = 2              B    
                           (6,4)    FC(row,col) =   0                R
                          :                             =>         :
	                      :                                        :
	               (257,3)          FC(row,col) = 2              B    
                           (258,4)  FC(row,col) =   0                R
                   
                   0,1,2 => R, G, B  for filters = 0x16161616 */ 
                for (c = FC(row,col); col < left+TS && col < width-3; col+=2) {
					pix = image + row*width+col;
					//val = ((pix[-1][FC(row,col-1)] + pix[0][c] + pix[1][FC(row,col+1)]) * 2 - pix[-2][c] - pix[2][c] + 2) >> 2; // horizontal
					val = (((pix[-1][1] + pix[0][c] + pix[1][1]) << 1) - pix[-2][c] - pix[2][c] + 2) >> 2; // horizontal
					/* pix[-1][1] = green channel of pixel on the left of current pixel
					      [-1] = position is current (row*width+col) minus 1 (-1)
					      [ 1] = color channel is 1 (Green)
					   >>2 => divided by 4   */
					val = CLIP(val);  // round val into [0,65535]
					rgb[0][row-top][col-left][1] = val;
					val = (((pix[-width][1] + pix[0][c] + pix[+width][1]) << 1) - pix[-2*width][c] - pix[2*width][c] + 2) >> 2; // vertical
					val = CLIP(val);
					rgb[1][row-top][col-left][1] = val;
					/* [0] = horizontal
					   [row-top][col-left] => img[3][3] to rgb[0][0], always in range[256][256]
					   [1] = channel 1 = channel Green */
				}
			}
			/*  Interpolate red and blue, and convert to CIELab:		*/
			for (d=0; d < 2; d++){ // dimension: 0 is horizontal and 1 is vertical for Green interpolation
				for (row=top+1; row < top+TS-1 && row < height-4; row++){       // row = 4 ~ 259, 260 ~ 515, ... 
					for (col=left+1; col < left+TS-1 && col < width-4; col++) { // col = 4 ~ 259, 260 ~ 515, ...
						pix = image + row*width+col;  // pix[position][color channel], raw data
						rix = &rgb[d][row-top][col-left]; // rgb, after Green interpolation
						lix = &lab[d][row-top][col-left]; // lab
						if ((c = 2 - FC(row,col)) == 1) { // if c = FC(row,col) = 1 = Green channel
							c = FC(row+1,col);// c = color pattern of the pixel below (row,col) where (row,col) is Green
							                  // so c is either Blue or Red
							val = pix[0][1] + (( pix[-1][2-c] + pix[1][2-c]- rix[-1][1] - rix[1][1] + 1) >> 1);  // horizontal (R or B)
							/*    G(0,0)          R(0,-1) , R(0,1) if c = 2   G^(0,-1)     G^(0,1)     
							                      B(0,-1) , B(0,1) if c = 0                        
							      if c = FC(row+1,col) = 2 , then val = R^(0,0) */
							// G^ => esitimated Green value
							val = CLIP(val);
							rix[0][2-c] = val; 							
							val = pix[0][1] + (( pix[-width][c] + pix[+width][c]- rix[-TS][1] - rix[+TS][1] + 1) >> 1);// vertical (B or R)
							/*    G(0,0)          B(-1,0) , B(1,0) if c = 2       G^(-1,0)      G^(1,0)     
							                      R(-1,0) , R(1,0) if c = 0                        
							      if c = FC(row+1,col) = 2 , then val = B^(0,0) */
					    }
	                    else { // if c = FC(row,col) = 0 or 2 => current point is Blue or Red
							   // in this situation, one channel of B and R is known, only need to compute another one.
							val = rix[0][1] + (( pix[-width-1][c] + pix[-width+1][c]+ pix[+width-1][c] + pix[+width+1][c]- rix[-TS-1][1] - rix[-TS+1][1] - rix[TS-1][1] - rix[TS+1][1] + 2) >> 2);
							//        G^(0,0)         	R(-1,-1) , R(-1,1) , R(1,-1) , R(1,1) if c = 2 - FC(row,col) = 0;       G^(-1,-1) , G^(-1,1) , G^(1,-1) , G^(1,1)       
					        //        if c = 0 , then FC(row,col) = 2, B(0,0) is known, the above equation is used to compute R^(0,0)
					        //        if c = 2 , then FC(row,col) = 0, R(0,0) is known, the above equation is used to compute B^(0,0)
							val = CLIP(val);
						}
						rix[0][c] = val;
						c = FC(row,col); 
						rix[0][c] = pix[0][c]; // fill one channel in (0,0) with known pixel value
						// For now, all pixels have been estimated in RGB in 2 directions, now start to translate to Lab
						// d = 0, Green is computed on horizontal direction, and R,B are estimated based on G^
						// d = 1, Green is computed on vertical direction
						
						// Then transfer into Lab 						
						xyz[0] = xyz[1] = xyz[2] = 0.5;
						FORC3 {
							xyz[0] += xyz_cam[0][c] * rix[0][c];
							xyz[1] += xyz_cam[1][c] * rix[0][c];
							xyz[2] += xyz_cam[2][c] * rix[0][c];
						}
						xyz[0] = cbrt[CLIP((int) xyz[0])];
						xyz[1] = cbrt[CLIP((int) xyz[1])];
						xyz[2] = cbrt[CLIP((int) xyz[2])];
						lix[0][0] = 64 * (116 * xyz[1] - 16);        // L
						lix[0][1] = 64 * 500 * (xyz[0] - xyz[1]);    // a
						lix[0][2] = 64 * 200 * (xyz[1] - xyz[2]);    // b
					}
				}
			}
			/*  Build homogeneity maps from the CIELab images: */
			memset (homo, 0, 2*TS*TS);
			for (row=top+2; row < top+TS-2 && row < height-5; row++) {
				tr = row-top;
				for (col=left+2; col < left+TS-2 && col < width-5; col++) {
					tc = col-left;
					// compute Lab difference
					for (d=0; d < 2; d++) {
						lix = &lab[d][tr][tc];
						for (i=0; i < 4; i++) { 
							l_diff[d][i] = ABS(lix[0][0]-lix[dir[i]][0]); // compute L difference on 4 directions
							ab_diff[d][i] = SQR(lix[0][1]-lix[dir[i]][1]) + SQR(lix[0][2]-lix[dir[i]][2]); // difference a & b
							// l_ab_diff[dimension][direction]
						}
					}
					// Adaptive Parameters (find the minimum possible boundary epsilon)
					l_eps = MIN(MAX(l_diff[0][0],l_diff[0][1]),  MAX(l_diff[1][2],l_diff[1][3]));
					ab_eps = MIN(MAX(ab_diff[0][0],ab_diff[0][1]),	MAX(ab_diff[1][2],ab_diff[1][3]));
					for (d=0; d < 2; d++){
						for (i=0; i < 4; i++){
							if (l_diff[d][i] <= l_eps && ab_diff[d][i] <= ab_eps){
								homo[d][tr][tc]++;   // higher homo => more similar => potential interpolate direction
							}
						}
					}
				}
			}
			/*  Combine the most homogenous pixels for the final result:	*/
			for (row=top+3; row < top+TS-3 && row < height-6; row++){
				tr = row-top;
				for (col=left+3; col < left+TS-3 && col < width-6; col++){
					tc = col-left;
					for (d=0; d < 2; d++){
						for (hm[d]=0, i=tr-1; i < tr+2; i++){
							for (j=tc-1; j < tc+2; j++) {
								hm[d] += homo[d][i][j];// compute homo direction in a 3x3 block
							}
						}
					}
					if (hm[0] != hm[1]) { // select direction with highest probability (homo value)
						FORC3 {image[row*width+col][c] = rgb[hm[1] > hm[0]][tr][tc][c];}
					}
					else{
						FORC3 {image[row*width+col][c] =(rgb[0][tr][tc][c] + rgb[1][tr][tc][c] + 1) >> 1;}
					/* fill the image with estimated value on the selected direction */
					}
				}
			}
		}
	}	
}																
								
int main()
{
	printf("Start.\n");
	FILE* fpin = NULL;
	printf("Open Binary Data File.\n");
	if ((fpin = fopen("data","r")) == NULL){
		fprintf(stderr,"fail.\n");
		return 0;
	}
	printf("Successfully Open Binary Data File.\n");
		
	raw_img = (ushort *) malloc (height*width*sizeof(ushort));	
	fread (raw_img,2,height*width,fpin); // read image data (16bit) as a 1D array
	printf("Successfully Load Data from Binary File.\n");	
	image = (ushort (*)[colors])calloc (height*width*sizeof(*image) , 2);	
	for (row=0;row<height;row++){
		for(col=0;col<width;col++){
			if (((row & 1) == 0) && ((col &1) == 0)){
				image[row*width+col][0] = raw_img[row*width+col];    // R (0,0)
			}
			if (((row & 1) == 0) && ((col &1) == 1)){
				image[row*width+col][1] = raw_img[row*width+col];;   // G1(0,1)
			}
			if (((row & 1) == 1) && ((col &1) == 1)){ 
				image[row*width+col][2] = raw_img[row*width+col];;   // B (1,1)
			}
			if (((row & 1) == 1) && ((col &1) == 0)){
				image[row*width+col][1] = raw_img[row*width+col];;   // G2(1,0)
			}
		}
	}
	printf("Successfully Transform to Image Shape.\n");
	printf("Shape of Image : %d x %d x %d\n",height,width,colors);
	printf("Simple Data Check : \n");
	for (row=0;row<1;row++){
		for(col=0;col<10;col++){
			for (color = 0;color<colors;color++){
				printf(" image[%d,%d][%d] = %d ",row,col,color,image[row*width+col][color]);
			}
			for (color = 0;color<colors;color++){
				
				if (image[row*width+col][color] == raw_img[row*width+col]){
					printf(" Same ");
				}
			}		
			printf("\n");
		}
	}			
	
	ahd_demosaic();
	printf("Demosaic Finished \n");
	printf("Start to Write Data.\n");
	FILE *write_ptr;
	write_ptr = fopen("Demosaic_Data","w+");  // w+ for write
	fwrite(image,2,width*height*colors,write_ptr); // num of points : width*height*colors
	                                               // each point has 2 bytes (16 bit)
	printf("Successfully Write Data.\n");
	fclose(fpin);
	fclose(write_ptr);
	free(raw_img);
	free(image);
	return 0;
}

