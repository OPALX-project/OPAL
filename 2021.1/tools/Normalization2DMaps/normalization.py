#!/usr/bin/python

import sys,re,os,string,fileinput

def printUsage():
    print "./normalization.py --fieldmap=FILE"

def nextLine(fh):
    line = fh.readline()
    if len(line) == 0:
        return line

    while line.startswith("#"):
        line = fh.readline()
    line = line.split("#")[0]
    return line.strip()

def computeRatio(fieldmap):
    fh = open(fieldmap, 'r')
    fileType = nextLine(fh)

    if not (fileType.startswith("2DMagnetoStatic") or
            fileType.startswith("2DElectroStatic") or
            fileType.startswith("2DDynamic")):
        print "./normalization.py only works for 2D field maps"
        sys.exit()


    orientation = fileType.split()[1].strip()
    fileType = fileType.split()[0].strip()

    dimension1 = nextLine(fh)
    if fileType.startswith("2DDynamic"):
        nextLine(fh)
    dimension2 = nextLine(fh)

    Nz = int(dimension1.split()[2]) + 1
    Nr = int(dimension2.split()[2]) + 1
    columnNr = 0

    if orientation == "ZX":
        Nr,Nz = Nz,Nr
        columnNr = 1


    maxGlobal = 0.0
    maxOnAxis = 0.0
    line = nextLine(fh)
    lineNr = 0
    while len(line) > 0:
        val = line.split()[columnNr].strip()
        val = abs(float(val))
        if val > maxGlobal:
            maxGlobal = val
        if ((orientation == "XZ" and lineNr < Nz) or
            (orientation == "ZX" and lineNr % Nr == 0)):
            if val > maxOnAxis:
                maxOnAxis = val
        line = nextLine(fh)
        lineNr += 1

    print "divide the amplitude by " + str(maxGlobal / maxOnAxis)

def main(argv):
    fieldmap = ""
    fieldmapSet = False
    for arg in argv:
        if arg.startswith("--fieldmap"):
            fieldmap = arg.split("=")[1]
            fieldmapSet = True
        elif arg.startswith("--help"):
            printUsage()
            sys.exit()

    if fieldmapSet:
        computeRatio(fieldmap)
    else:
        printUsage()

#call main
if __name__ == "__main__":
    main(sys.argv[1:])