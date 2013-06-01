#!/usr/bin/env python
import argparse
import numpy as np
import os
from subprocess import call

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def clean_files(inp_dict):
    os.remove(inp_dict["inputname"])
    os.remove(inp_dict["xyzw"])

def combine_files(inp_dict, start, stop, nPoints):
    nxyz = file_len(inp_dict["xyzw"]) 
    print nxyz
    fxyzw = open(inp_dict["xyzw"], 'r')
    grid = np.linspace(start, stop, nPoints)
    f = open(inp_dict["gridfile"], 'w')
    f.write("{ntot:10d}{nint:10d}\n".format(ntot=nxyz+nPoints, nint=nxyz+1))
    for line in fxyzw.readlines():
        f.write(line)
    fxyzw.close()
    for i in range(len(grid)):
        f.write("{x:25.14E}{y:25.14E}{z:25.14E}{w:25.14E}\n".format(x=0.0, y=0.0, z=grid[i], w=0.0))
    f.close()

def get_input_variables(anglo, basinfo, debug, gridfile, nrad, pruning):

    inp_dict = {}
    inp_dict["inputname"]  = "grid.input"
    inp_dict["xyzw"]       = "xyzw.xyzw"
    inp_dict["anglo"]      = str(anglo)
    inp_dict["basinfo"]    = basinfo
    if debug:
        inp_dict["debug"]  = ".true."
    else:
        inp_dict["debug"]  = ".false."
    inp_dict["gridfile"]   = gridfile
    inp_dict["nradial"]    = str(nrad)
    if pruning:
        inp_dict["pruning"] = ".true."
    else:
        inp_dict["pruning"] = ".false."
    return inp_dict

def write_grid_input(inp):
    '''Write the input file for the grid program.'''

    input = open(inp["inputname"], 'w')
    input.write("&input\n")
    input.write("angularLOrder={0:s}\n".format(inp["anglo"]))
    input.write("basinfo='{0:s}'\n".format(inp["basinfo"]))
    input.write("debug={0:s}\n".format(inp["debug"]))
    input.write("gridfile='{0:s}'\n".format(inp["xyzw"]))
    input.write("nRadialPoints={0:s}\n".format(inp["nradial"]))
    input.write("pruning={0:s}\n".format(inp["pruning"]))
    input.write("/")
    input.close()

def get_integration_grid(inp):
    write_grid_input(inp) 
    grid_program = '/home/lmentel/work/mercurial_repo/Vconstr/Source/generate_grid/Grid.x'
    p = call([grid_program, inp["inputname"]])
    print p

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("basinfo",
                        help = "basinfo file")
    parser.add_argument("gridfile",
                        help = "output grid file")
    parser.add_argument("start",
                        type=float,
                        help="start")
    parser.add_argument("stop",
                        type=float,
                        help="stop")
    parser.add_argument("nPoints",
                        type=int,
                        help="number of points")

    parser.add_argument("-a",
                        "--anglo",
                        type=int,
                        choices = [5,7,11,13,17,19,21,23,29,31,35,41,47,53,59], 
                        default = 11,
                        help = "angular L Order")
    parser.add_argument("-nrad",
                        type = int,
                        default = 50, 
                        help = "number of radial points")
    parser.add_argument("-p",
                        "--pruning",
                        action = "store_true",
                        help = "pruning")
    parser.add_argument("-d",
                        "--debug",
                        action = "store_true",
                        help = "debug")
    args = parser.parse_args()

    inp_dict = get_input_variables(args.anglo, args.basinfo, args.debug, args.gridfile, 
                                   args.nrad, args.pruning)

    get_integration_grid(inp_dict)

    combine_files(inp_dict, args.start, args.stop, args.nPoints)

    clean_files(inp_dict)

if __name__ =="__main__":
    main()
