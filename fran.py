#! /usr/bin/env python
# coding: utf-8

## @file fran.py
## @author Cristian A. Vega Martinez <cnvega(at)fcaglp.unlp.edu.ar>
## @copyright Copyright © 2017 Cristian A. Vega M.
##
## @brief Fractal class
## 
## FrAn - Fractal Analyzer software.
## The fractal class allows to load an RBG image and can calculate its fractal spectrum, 
## the null model (randomizing the image) and the corresponding fractal excess.

import sys
import os
from PIL import Image
from pyx import * # plots
#import psyco
import franTools
import imrand
#from readcol import *

#psyco.full()
text.set(lfs="foils17pt")

class Fractal:
    """A Fractal Image"""
    def __init__(self):
        self.filename = ''
        self.datafilename = ''
        self.origImage = None
        self.grayImage = None
        self.nullImage = None
        self.average = 0.0
        self.desvest = 0.0
        self.spectrum = None
        self.nullSpectrum = None
        self.excess = None
    
    def __set_atr(self):
        if self.origImage != None:
            self.grayImage = self.origImage.convert("L")
            self.average = franTools.pixelAverage(self.grayImage)
            self.desvest = int(franTools.pixelDesvest(self.grayImage))

    def load(self, filename):
        try: 
            self.origImage = Image.open(filename)
        except IOError:
            print("Cannot load image file")
            return
        self.filename = filename
        self.__set_atr()
    
    def loadImage(self, imagein):
        self.origImage = imagein
        self.__set_atr()

    def fractalDimension(self, thr = -1):
        if thr == -1:
            thr = self.average
        result, data = franTools.fracDim(self.grayImage, thr)
        if result == 0:
            print("Data out of range!")
            return 0, 0
        if result == -1:
            print("Division by zero!")
            return 0, 0
        for element in result:
            print("D_"+element[0]+" = "+str(element[1])+"\t\terr = "+str(element[3]))
        return result[0][1], result[0][3]

    def fractalSpectra(self, specPoints = 31):
        dev = list(range(specPoints))
        thresVect = list(range(specPoints))
        spec = []
        specerr = []
        rango = 1.5
        # Checking range limits:
        if self.average + 1.5*self.desvest > 255:
            rango = float(255 - self.average)/float(self.desvest)
        if self.average - 1.5*self.desvest < 0:
            rango = float(self.average)/float(self.desvest)
        # Processing fractal dimension at different points:
        for i in range(specPoints):
            dev[i] = -rango+i*2*rango/(specPoints-1.0)
            thresVect[i] = int(round(self.average + self.desvest*dev[i]))
            print("--------------------------------------")
            print("%%% Fractal Dimension with threshold: "+str(thresVect[i]))
            auxFD, auxFDerr = self.fractalDimension(thresVect[i])
            if auxFD != 0:
                spec.append(auxFD)
                specerr.append(auxFDerr)
            else:
                spec.append(0.)
                specerr.append(0.)
            #print(spec) 
        self.spectrum = thresVect, dev, spec, specerr
        return thresVect, dev, spec, specerr
    
    def nullModel(self, randStep = 20):
        # Randomization constants:
        condition = 0.001
        imMaxX, imMaxY = self.grayImage.size
        maxRandPix = imMaxX+imMaxY
        loops = maxRandPix*randStep/100.
        if (loops-int(loops)) > 0.: # rounding
            loops = int(loops) + 1
        else: 
            loops = int(loops)
        # First randomization - 100%:
        randImage = self.grayImage.copy()
        print("Randomizing the image...")
        for i in range(int(maxRandPix)):
            imrand.randomCol(randImage)
            imrand.randomRow(randImage)
        aux = Fractal()
        aux.loadImage(randImage.copy())
        result = aux.fractalSpectra(10)
        oldSp = result[2]
        del aux, result
        # After 100% spectra variability test:
        different = True
        while different == True:
            print("Randomizing...")
            for i in range(loops):
                imrand.randomCol(randImage)
                imrand.randomRow(randImage)
            nueva = Fractal()
            nueva.loadImage(randImage.copy())
            result = nueva.fractalSpectra(10)
            newSp = result[2]
            meanSp = 0.0
            for i in range(len(newSp)):
                meanSp += (newSp[i] - oldSp[i])
            meanSp = abs(meanSp)/10.
            if meanSp <= condition:
                different = False
                self.nullImage = randImage
                self.nullSpectrum = nueva.fractalSpectra()
            else:
                oldSp = newSp
            del nueva, result
        return self.nullSpectrum

    def excessSpectra(self):
        if self.spectrum == None:
            self.fractalSpectra()

        if self.nullSpectrum == None:
            self.nullModel()
        
        self.excess = []
        for i in range(31):
            aux = self.nullSpectrum[2][i] - self.spectrum[2][i]
            self.excess.append(aux)

        return self.excess

    def saveData(self):
        if self.excess == None:
            self.excessSpectra()

        mainName = os.path.splitext(self.filename)[0]
        datafile = open(mainName+".dat","w") 
        datafile.write("# thres\tdev\tspec\tspec_err\tnullsp\texcess\n")
        for i in range(31):
            for j in range(4):
                datafile.write(str(self.spectrum[j][i])+"\t")
            datafile.write(str(self.nullSpectrum[2][i])+"\t")
            datafile.write(str(self.excess[i])+"\n")

        self.grayImage.save(mainName+"_gray.png","PNG")
        self.nullImage.save(mainName+"_null.png","PNG")
        self.datafilename = mainName+".dat"
    
    def graphExcess(self):
        if self.datafilename == '':
            self.excessSpectra()
            self.saveData()
        # Fractal excess spectra graphic:
        mainName = os.path.splitext(self.filename)[0]
        graph1 = graph.graphxy(width=16,
            x=graph.axis.linear(title="Relative threshold value"),
            y=graph.axis.linear(title="Fractal Excess"))
        graph1.plot(graph.data.file(self.datafilename,x=2,y=6),
            [graph.style.symbol(symbol=graph.style.symbol.triangle,size=0.3,
            symbolattrs=[deco.filled([color.rgb.blue])])])
        graph1.writeEPSfile(mainName+"_excess")
        # Fractal Dimension and null model graphics:
        graph2 = graph.graphxy(width=16,
            x=graph.axis.linear(title="Relative threshold value"),
            y=graph.axis.linear(title="Fractal Dimension"),
            key=graph.key.key(pos="br",dist=0.1))
        graph2.plot(graph.data.file(self.datafilename,x=2,y=3,title="original"),
            [graph.style.symbol(symbol=graph.style.symbol.square,size=0.3,
            symbolattrs=[deco.filled([color.rgb.green])])])
        graph2.plot(graph.data.file(self.datafilename,x=2,y=5,title="randomized"),
            [graph.style.symbol(symbol=graph.style.symbol.diamond,size=0.3,
            symbolattrs=[deco.filled([color.rgb.red])])])
        graph2.writeEPSfile(mainName+"_specs")

if __name__ == "__main__":
    myfractal = Fractal()
    myfractal.load(sys.argv[1])
    myfractal.excessSpectra()
    myfractal.saveData()
    myfractal.graphExcess()

