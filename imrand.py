#! /usr/bin/env python
# coding: utf-8

## @file imrand.py
## @author Cristian A. Vega Martinez <cnvega(at)fcaglp.unlp.edu.ar>
## @copyright Copyright © 2017 Cristian A. Vega M.
##
## @brief Null model: randomize images.
## 
## FrAn - Fractal Analyzer software.
## The functions defined here are useful to find the null models of the images
## from which the fractal excess are going to be calculated.

import os, sys
from PIL import Image
from numpy import random

def randomPix(im):
    """Switch the value of two random pixels"""
    x1 = random.random_integers(0,im.size[0]-1)
    x2 = random.random_integers(0,im.size[0]-1)
    y1 = random.random_integers(0,im.size[1]-1)
    y2 = random.random_integers(0,im.size[1]-1)

    aux = im.getpixel((x1,y1))
    im.putpixel((x1,y1),im.getpixel((x2,y2)))
    im.putpixel((x2,y2),aux)
    return 1


def randomCol(im):
    """Switch two random columns in the image"""
    x1 = random.random_integers(0,im.size[0]-1)
    x2 = random.random_integers(0,im.size[0]-1)
    
    col1 = (x1,0,x1+1,im.size[1])
    col2 = (x2,0,x2+1,im.size[1])
    aux = im.crop(col1)
    aux.load()
    im.paste(im.crop(col2),col1)
    im.paste(aux,col2)
    return 1


def randomRow(im):
    """Switch two random rows in the image"""
    y1 = random.random_integers(0,im.size[1]-1)
    y2 = random.random_integers(0,im.size[1]-1)
    
    row1 = (0,y1,im.size[0],y1+1)
    row2 = (0,y2,im.size[0],y2+1)
    aux = im.crop(row1)
    aux.load()
    im.paste(im.crop(row2),row1)
    im.paste(aux,row2)
    return 1


def randomizar(option, infile, maxPercent, delta):
    try:    # Loading image file
              imName = Image.open(infile)
    except IOError:
              print("File doesn't exist!")
              return 0
    
    fileName = os.path.splitext(infile)[0]
    imGray = imName.convert("L")
    imMaxX, imMaxY = imName.size # Max sizes of the image
    maxPercent = int(maxPercent) # converts string to int 
    delta = float(delta)

    if option == '-p':
        maxRandPix = imMaxX*imMaxY   # 100% 
    elif option == '-c':
        maxRandPix = imMaxX+imMaxY   # 100%

    loops = maxRandPix*maxPercent/100.  # pixel switching number
    deltaPix = int(maxRandPix*delta/100.)    # delta switching
    j=0.

    if (loops-int(loops)) > 0.:
        loops = int(loops) + 1
    else: 
        loops = int(loops)

    print(("deltaPix = "+str(deltaPix)))
    if option == '-p':
        for i in range(loops):
            totalPer = int((i+1.)*100./loops)
            percentPix = int((i+1.)*100./maxRandPix)
            sys.stdout.write("  Processing... "+str(totalPer)+"%\r")
            sys.stdout.flush() 
            randomPix(imGray)
            if (i+1.) >= deltaPix*(j+1.):  # Saving partial images
                imGray.save(fileName+"_"+str(percentPix)+"prand.png","PNG")
                j+=1.

    if option == '-c':
        for i in range(loops):
            totalPer = int((i+1.)*100./loops)
            percentPix = int((i+1.)*100./maxRandPix)+1 # why is this +1 here?
            sys.stdout.write("  Processing... "+str(totalPer)+"%\r")
            sys.stdout.flush()
            randomCol(imGray)
            randomRow(imGray)
            if (i+1.) >= deltaPix*(j+1.):  # Saving partial images
                imGray.save(fileName+"_"+str(percentPix)+"crrand.png","PNG")
                j+=1.
    
    print("Done!")
    return 1

if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == '-p':
        randomizar(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    elif len(sys.argv) == 5 and sys.argv[1] == '-c':
        randomizar(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("""Usage: imrand OPTION FILENAME MR DR
ImRand - Image randomizer
  MR    Max percent for randomize
  DR    Delta for save each image
OPTIONS:
  -c    randomize using column-row mode
  -p    randomize using pixel-pixel mode
        """)

