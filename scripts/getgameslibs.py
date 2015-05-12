#!/usr/bin/env python

import argparse
import os
import re
import shutil
import subprocess

def main():
    '''
    Create the libgamess.a and copy libgamess.a and libddi.a to a new dir
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gmspath",
                        help="Path to GAMESS(US) installation (GMSPATH)")
    parser.add_argument("-d", "--destination",
                        default=os.getcwd(),
                        help="Destination path where the lib will be copied")
    args = parser.parse_args()

    destination = os.path.abspath(args.destination)

    if args.gmspath:
        gmspath = os.path.abspath(args.gmspath)
    else:
        if os.getenv("GMS_PATH") is not None:
            gmspath = os.getenv("GMS_PATH")
        elif os.getenv("GMSPATH") is not None:
            gmspath = os.getenv("GMSPATH")
        raise ValueError("gmspath needs to be specified")

    comp = os.path.join(gmspath, "comp")
    objdir = os.path.join(gmspath, "object")
    libddi = os.path.join(gmspath, "ddi", "libddi.a")
    gmssrc = os.path.join(gmspath, "source", "gamess.src")
    libgamess = os.path.join(objdir, "libgamess.a")

    with open(gmssrc, 'r') as gms:
        src = gms.read()

    src = re.sub(r'\s{7}PROGRAM GAMESS', '\n      SUBROUTINE GAMESS', src)
    src = re.sub(r'\s{7}END\n', '\n      END SUBROUTINE\n', src, count=1)

    # rename the gamess.src file
    gmsorig = os.path.join(gmspath, "source", "gamess.src.orig")
    os.rename(gmssrc, gmsorig)

    with open(gmssrc, 'w') as gms:
        gms.write(src)

    # compile the new gamess.src file
    os.chdir(gmspath)
    subprocess.call([comp, "gamess"])

    # go to the objects dir and create libgamess.a
    os.chdir(objdir)
    subprocess.call("ar crv libgamess.a *.o", shell=True)

    # switch back to the original gamess.src file
    os.rename(gmsorig, gmssrc)

    shutil.copyfile(libgamess, os.path.join(destination, "libgamess.a"))
    shutil.copyfile(libddi, os.path.join(destination, "libddi.a"))

if __name__ == "__main__":
    main()
