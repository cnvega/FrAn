#! /usr/bin/env python

import sys

def readcol(filename, colnumber):
    """Reads a column from a file"""
    vector = []
    try:
        thefile = open(filename,"r")
        for linea in thefile.readlines():
            if linea[0] != '#':
                vector.append(linea.split()[int(colnumber)-1])
    except:
        return 0
    return vector

def readcols(filename, *colnumbers):
    """Reads columns from a file"""
    vector = []
    try:
        for numbers in colnumbers:
            vector.append([])
        thefile = open(filename,"r")
        for linea in thefile.readlines():
            if linea[0] != '#':
                for aux in range(len(colnumbers)):
                    try:
                        vector[aux].append(linea.split()[int(colnumbers[aux])-1])
                    except:
                        pass
    except:
        return 0
    return vector

if __name__ == "__main__":
    if len(sys.argv) == 3:
        vector = readcol(*sys.argv[1:])
        if vector != 0:
            for element in vector:
                print(element)
    else:
        print("""Usage: readcol N""")

