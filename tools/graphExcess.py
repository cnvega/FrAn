#! /usr/bin/env python

import sys
import os
from pyx import * # plots

text.set(lfs="foils17pt")
def graphExcess(datafilename):
    # Fractal excess spectra graphic:
    mainName = os.path.splitext(datafilename)[0]
    graph1 = graph.graphxy(width=16,
        x=graph.axis.linear(title="Relative thresholding value"),
        y=graph.axis.linear(title="Fractal Excess"))
    graph1.plot(graph.data.file(datafilename,x=2,y=6),
        [graph.style.symbol(symbol=graph.style.symbol.triangle,size=0.3,
        symbolattrs=[deco.filled([color.rgb.blue])])])
    graph1.writeEPSfile(mainName+"_excess")
    # Fractal Dimension and null model graphics:
    graph2 = graph.graphxy(width=16,
        x=graph.axis.linear(title="Relative thresholding value"),
        y=graph.axis.linear(title="Fractal Dimension"),
        key=graph.key.key(pos="br",dist=0.1))
    graph2.plot(graph.data.file(datafilename,x=2,y=3,title="original"),
        [graph.style.symbol(symbol=graph.style.symbol.square,size=0.3,
        symbolattrs=[deco.filled([color.rgb.green])])])
    graph2.plot(graph.data.file(datafilename,x=2,y=5,title="randomized"),
        [graph.style.symbol(symbol=graph.style.symbol.diamond,size=0.3,
        symbolattrs=[deco.filled([color.rgb.red])])])
    graph2.writeEPSfile(mainName+"_specs")

if __name__ == "__main__":
    graphExcess(sys.argv[1])
