#! /usr/bin/env python
# coding: utf-8

## @file franTools.py
## @author Cristian A. Vega Martinez <cnvega(at)fcaglp.unlp.edu.ar>
## @copyright Copyright © 2017 Cristian A. Vega M.
##
## @brief Fractal Analyzer Toolkit
##
##  This is the toolkit of FrAn - Fractal Analyzer software.
## These functions can take an image independently of the file type and apply 
## a filter to convert it from RGB to greyscale according to a given threshold. 
## Then, using the box counting method, the Fractal Dimension or the Fractal 
## Spectrum can be obtained depending on the given option.

import sys
from PIL import Image
from math import log, exp
from pyx import * # Module to make plots
#import psyco
import copy

#psyco.full()


def thresholding(imageIn,cut): 
    """Applies threshold filter to 2D gray image using cut value."""
    result = imageIn.copy()
    result = result.point(lambda i: 255)
    zeros = result.point(lambda i: i*0)
    mask = imageIn.point(lambda i: i < cut and 255)
    result.paste(zeros, None, mask)
    return result

def boxCounting2(imageIn, boxSize):
    whites = 0.0
    blacks = 0.0
    both = 0.0
    nx = int(imageIn.size[0]/boxSize)
    ny = int(imageIn.size[1]/boxSize)
    boxlist = []
    addbox = boxlist.append
    coords = [0,0,0,0,0,0,0,0]
    x0,y0,x1,y1,x2,y2,x3,y3 = 0,1,2,3,4,5,6,7    
    for ix in range(nx):
        for iy in range(ny):
            coords[x0] = coords[x1] = boxSize*ix
            coords[y0] = coords[y3] = boxSize*iy
            coords[y1] = coords[y2] = boxSize*iy + boxSize
            coords[x2] = coords[x3] = boxSize*ix + boxSize
            boxlist.append(copy.copy(coords))
    for box in boxlist:
        boxImage = imageIn.transform((boxSize,boxSize), Image.QUAD, box)
        histogram = boxImage.histogram()
        if histogram[0] == 0:
            whites += 1.0
        elif histogram[255] == 0:
            blacks += 1.0
        else:
            both += 1.0
    return whites, blacks, both, nx*ny


def boxCounting(imageIn, boxSize): 
    """ Returns a tuple with the results of the box counting method using
    boxSize as size of each box. (whites, blacks, both, N)"""
    # Box shape:
    #  (x0,y0)----(x3,y3)
    #     |          |
    #     |          |
    #  (x1,y1)----(x2,y2)
    whites = 0.0
    blacks = 0.0
    both = 0.0
    nx = int(imageIn.size[0]/boxSize)
    ny = int(imageIn.size[1]/boxSize)
    coords = [0,0,0,0,0,0,0,0]
    x0,y0,x1,y1,x2,y2,x3,y3 = 0,1,2,3,4,5,6,7    
    for ix in range(nx):
        for iy in range(ny):
            coords[x0] = coords[x1] = boxSize*ix
            coords[y0] = coords[y3] = boxSize*iy
            coords[y1] = coords[y2] = boxSize*iy + boxSize
            coords[x2] = coords[x3] = boxSize*ix + boxSize
            boxImage = imageIn.transform((boxSize,boxSize), Image.QUAD, coords)
            histogram = boxImage.histogram()
            if histogram[0] == 0:
                whites += 1.0
            elif histogram[255] == 0:
                blacks += 1.0
            else:
                both += 1.0

    return whites, blacks, both, nx*ny


def pixelAverage(imageIn): # Calculates an average of a greyscale image
    """Returns an average of all the pixels in the image"""
    histogram = imageIn.histogram()
    suma = 0.0
    result = 0.0
    for counter in range(256):
        suma += histogram[counter]
        result += histogram[counter]*counter

    return int(result/suma)


def pixelDesvest(imageIn): # Calculates the standard deviation of an image
    """Returns the standard desviation of all the pixels in the image"""
    histogram = imageIn.histogram()
    mean = pixelAverage(imageIn)
    suma = 0.0
    total = 0.0
    for counter in range(256):
        total += histogram[counter]
        suma += histogram[counter]*((counter-mean)**2)

    return (suma/total)**0.5


def linFit(xvec, yvec): # Make a linear fit of the points (x,y)
    """Returns a linear fit of the X and Y vector using Y = mX + n.
    This function returns m, n, dm and dn in a tuple."""
    sumx = 0.0
    sumy = 0.0
    sumxx = 0.0
    sumxy = 0.0
    total = len(xvec)
    for i in range(len(xvec)):
        sumx += xvec[i]
        sumy += yvec[i]
        sumxx += (xvec[i])**2.0
        sumxy += (xvec[i])*(yvec[i])

    delta = total*sumxx - sumx**2.0
    m = (total*sumxy-sumx*sumy)/delta
    n = (sumxx*sumy - sumx*sumxy)/delta
    aux = 0.0
    for i in range(len(xvec)):
        aux += (yvec[i] - (m*xvec[i] + n))**2.0
    
    sigma = aux/(total-2.0)
    dm = (sigma*total/delta)**0.5
    dn = (sigma*sumxx/delta)**0.5

    return m,n,dm,dn    


def unique(thelist):
    """Returns the list without repeated elements"""
    return list(set(thelist))


def imageProp(imageName):
    """Shows image properties"""
    originalImage = Image.open(imageName)
    print("Original Image Properties:")
    print("Image name: "+str(imageName))
    print("Image type: "+str(originalImage.mode))
    print("Image size: "+str(originalImage.size))
    grayImage = originalImage.convert("L")
    average = pixelAverage(grayImage)
    desvest = int(pixelDesvest(grayImage))
    print("Grayscale Image Stats:")
    print("  mean:  "+str(average)+"\n  sigma: "+str(desvest))
    return 0


def fracDim(imageIn, thres):
    """Obtains the fractal dimension at a certain thresholding."""
    print("Processing Fractal Dimension...")
    if imageIn.mode != "L":
        print("Applying grayscale filter...")
        imageIn = imageIn.convert("L")

    print("Applying Thresholding ["+str(thres)+"]")
    bwImage = thresholding(imageIn, thres)
    wholeHistogram = bwImage.histogram()
    if wholeHistogram[0] == 0 or wholeHistogram[255] == 0:
        return 0, 0
    numSteps = 30
    area = (bwImage.size[0]*bwImage.size[1])**0.5
    minStep = int(0.001*area)  # 0.1%
    maxStep = int(0.3*area)    # 30%
    if minStep < 2:
        minStep = 2
    
    print("Generating size boxes vector")
    allData, result = [], []
    lnr, lnBW, lnWBW, lnBBW = [], [], [], []
    lnMinStep = log(1.0/minStep)
    lnMaxStep = log(1.0/maxStep)
    lnStepSize = (lnMaxStep-lnMinStep)/(numSteps-1.0)
    boxesVect = list(range(numSteps))
    for i in range(numSteps):
        boxesVect[i] = int(round(1.0/exp(lnMinStep+i*lnStepSize)))
    
    boxesVect = unique(boxesVect)
    boxesVect.sort()
    print(boxesVect)
    for box in boxesVect:
        W,B,BW,N =  boxCounting(bwImage, box)
        try:
            logBW = log(BW)
            logWBW = log(W+BW)
            logBBW = log(B+BW)
            lnBW.append(logBW)
            lnWBW.append(logWBW)
            lnBBW.append(logBBW)
            lnr.append(log(1.0/box))
            allData.append((W,B,BW,N,box,log(1.0/box),logBW,logWBW,logBBW))
        except OverflowError:
            print("Infinity")

    try:
        aux = linFit(lnr, lnBW)
        result.append(("BW",aux[0],aux[1],aux[2],aux[3]))
        aux = linFit(lnr, lnWBW)
        result.append(("WBW",aux[0],aux[1],aux[2],aux[3]))
        aux = linFit(lnr, lnBBW)
        result.append(("BBW",aux[0],aux[1],aux[2],aux[3]))
    except ZeroDivisionError:
        return -1, 0

    return result, allData


def fran(imageName, option): ### Main of the script

    print("FrAn - Fractal Analyzer Software")
    print("Loading image...")
    originalImage = Image.open(imageName) # image load
    print("Image type: "+str(originalImage.mode))
    print("Image size: "+str(originalImage.size))
    print("Applying greyscale filter...")
    grayImage = originalImage.convert("L")
    average = pixelAverage(grayImage)
    desvest = int(pixelDesvest(grayImage))
    print("Image Stats:")
    print("  Mean:  "+str(average)+"\n  sigma: "+str(desvest))
    text.set(lfs="foils17pt")
    # Using -fd or -t option (fractal dimension at fixed threshold):
    if option != "fs":
        print("% Fractal Dimension Analysis %")
        if option == "fd":
            thr = average
            outName = imageName.replace(".","")+"_FD"  # output file
        else:
            thr = int(option)
            outName = imageName.replace(".","")+"_Ft"+option  # output file
        fractal, data = fracDim(grayImage, thr)   # Fractal analysis
        if fractal == 0:
            print("Data out of range!")
            return 0
        if fractal == -1:
            print("Division by zero!")
            return 0
        thresImage = thresholding(grayImage, thr) # Saving thres. image
        thresImage.save(imageName+"_t"+str(thr)+".png", "PNG")
        datafile = open(outName+".dat","w")
        # file header: results
        datafile.write("# Image Name   = "+str(imageName)+"\n") 
        datafile.write("# Threshold    = "+str(thr)+"\n#\n")
        for element in fractal:
            datafile.write("# D_"+element[0]+" = "+str(element[1]))
            datafile.write("\t\terr = "+str(element[3])+"\n")
            datafile.write("# k_"+element[0]+" = "+str(element[2]))
            datafile.write("\t\terr = "+str(element[4])+"\n#\n")
        # file body: data
        datafile.write("#W\tB\tBW\tTotal\tr\tln(1/r)\tln(BW)\tln(W+BW)\tln(B+BW)\n")
        for element in data:
            for i in range(9):
               datafile.write(str(element[i])+"\t")
            datafile.write("\n")
        # closing file
        datafile.close()
        # graphic: Data and linear fits in an EPS file.    
        graphic=graph.graphxy(width=16,x=graph.axis.linear(title="ln(1/r)"),
            y=graph.axis.linear(title="ln(N)"), key=graph.key.key(pos="br",dist=0.1))
        graphic.plot([graph.data.file(outName+".dat", x=6, y=7, title="BW data"),
            graph.data.file(outName+".dat", x=6, y=8, title="W+BW data"),
           graph.data.file(outName+".dat", x=6, y=9, title="B+BW data")],
            [graph.style.symbol(size=0.3)])
        graphic.plot(graph.data.function("y(x)="+str(fractal[0][1])+"*x+"+str(fractal[0][2]),title="BW fit"),
            [graph.style.line([style.linestyle.solid, color.rgb.green])])
        graphic.plot(graph.data.function("y(x)="+str(fractal[1][1])+"*x+"+str(fractal[1][2]),title="W+BW fit"),
            [graph.style.line([style.linestyle.solid, color.rgb.red])])
        graphic.plot(graph.data.function("y(x)="+str(fractal[2][1])+"*x+"+str(fractal[2][2]),title="B+BW fit"),
            [graph.style.line([style.linestyle.solid, color.rgb.blue])])
        graphic.writeEPSfile(outName)
        graphic.writePDFfile(outName)
    # Using -fs option (fractal spectrum around average):
    if option == "fs":
        print("% Fractal Spectrum Analysis %")
        specPoints = 31
        dev = list(range(specPoints))
        thresVect = list(range(specPoints))
        outName = imageName.replace(".","")+"_FSp"    # output file
        datafile = open(outName+".dat","w")
        spectrum = []
        rango = 1.5
        if average+1.5*desvest > 255:
            rango = float(255-average)/float(desvest)
        if average-1.5*desvest < 0:
            rango = float(average)/float(desvest)
        for i in range(specPoints):
            dev[i] = -rango+i*2*rango/(specPoints-1.0)
            thresVect[i] = int(round(average + desvest*dev[i]))
            print("--------------------------------------")
            print("%%% Fractal Dimension with threshold: "+str(thresVect[i]))
            fractal, data = fracDim(grayImage, thresVect[i])  # Fractal analysis 
            if fractal != 0 and fractal != -1:
                aux = (thresVect[i], dev[i])
                spectrum.append(aux+fractal[0][1:]+fractal[1][1:]+fractal[2][1:])
        # Saving max/min threshold images
        thresImage = thresholding(grayImage, thresVect[0])
        thresImage.save(imageName+"_t"+str(thresVect[0])+".png", "PNG")
        thresImage = thresholding(grayImage, average)
        thresImage.save(imageName+"_t"+str(average)+".png", "PNG")
        thresImage = thresholding(grayImage, thresVect[30])
        thresImage.save(imageName+"_t"+str(thresVect[30])+".png", "PNG")
        # File header: name of each column
        datafile.write("# % Fractal Spectrum Analysis %\n")
        datafile.write("# Image Name: "+str(imageName)+"\n")
        datafile.write("# Thres\tdispl\t")
        datafile.write("D_bw\tk_bw\tD_err\tk_err\t")
        datafile.write("D_wbw\tk_wbw\tD_err\tk_err\t")
        datafile.write("D_bbw\tk_bbw\tD_err\tk_err\n")
        # File body: Fractal Spectrum data
        for element in spectrum:
            for counter in range(14):
                datafile.write(str(element[counter])+"\t")
            datafile.write("\n")
        # closing file
        datafile.close()
        # Graphic: BW, W-BW and B-BW fractal spectrums in the same graphic.
        graphic = graph.graphxy(width=16, 
            x=graph.axis.linear(title="Relative threshold value"),
            y=graph.axis.linear(title="Fractal Dimension"),
            key=graph.key.key(pos="br",dist=0.1))
        graphic.plot(graph.data.file(outName+".dat",x=2,y=3,dy=5,title="BW"),
            [graph.style.symbol(symbol=graph.style.symbol.square, size=0.3, 
            symbolattrs=[deco.filled([color.rgb.green])]), 
            graph.style.errorbar(size=0, errorbarattrs=[style.linewidth.thin])])
        graphic.plot(graph.data.file(outName+".dat", x=2, y=7, dy=9, title="W+BW"),
            [graph.style.symbol(symbol=graph.style.symbol.triangle, size=0.24, 
            symbolattrs=[deco.filled([color.rgb.red])]), 
            graph.style.errorbar(size=0, errorbarattrs=[style.linewidth.thin])])
        graphic.plot(graph.data.file(outName+".dat", x=2, y=11, dy=13, title="B+BW"),
            [graph.style.symbol(symbol=graph.style.symbol.diamond, size=0.24, 
            symbolattrs=[deco.filled([color.rgb.blue])]), 
            graph.style.errorbar(size=0, errorbarattrs=[style.linewidth.thin])])
        graphic.writeEPSfile(outName)
        graphic.writePDFfile(outName)

    print("Done!")
    
if __name__ == "__main__":
    if len(sys.argv) != 1:
        option = sys.argv[1] 
    else: 
        option = ''
    if option == '-fd' or option == '-d' or option == 'fd':
        fran(sys.argv[2],"fd")
    elif option == '-fs' or option == '-s' or option == 'fs':
        fran(sys.argv[2], "fs")
    elif option == '-t':
        fran(sys.argv[3],sys.argv[2])
    elif option == '-p':
        imageProp(sys.argv[2])
    else:
        print("""Usage: franTools OPTION InputImage
FrAn - Fractal Analyzer Software (Toolkit)
OPTIONS:
-fd,  -d    Calculates the fractal dimension using box counting
 fd         method and threshold at average pixel value.

-fs,  -s    Calculates the fractal spectrum with thresholds
 fs         around the average pixel value in 1.5 desvest units.

-t N        Calculates the fractal dimension using box counting
            method with threshold at N, where 0 < N < 255.

-p          Shows image properties only.
        """)

