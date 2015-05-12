#!/usr/bin/env python

import argparse
import os
import re
import sys

def get_info(lines, rege):
    '''Get the infromation form contents based on the regular expression.'''

    cp = re.compile(rege)

    for line in lines:
        match = cp.search(line)
        if match:
            return int(match.group(1))

def write_input(input_name, input_dict):
    inp = open(input_name, 'w')
    # write the input namelist
    inp.write("&input\n")
    for key in sorted(input_dict.iterkeys()):
        inp.write("\t{:s}={:s}\n".format(key, str(input_dict[key])))
    inp.write("/\n")
    # write the occupations namelist
    inp.write("&occupations\n")
    inp.write("\t{:s} = {:s}\n".format('occmo(1)', ", ".join(['2.0']*input_dict['nmos'])))
    inp.write("/\n")
    inp.close()

def get_input_data(logfile):
    '''Retrieve the available information form gamess-us logfile and put
    all the DMFT input keywords in a dictionary with some default values
    set.'''

    log = open(logfile, 'r')
    contents = log.readlines()
    log.close()

    re_nocc_alpha = r'NUMBER OF OCCUPIED ORBITALS \(ALPHA\)\s+=\s*(\d+)'
    re_nocc_beta  = r'NUMBER OF OCCUPIED ORBITALS \(BETA \)\s+=\s*(\d+)'

    na = get_info(contents, re_nocc_alpha)
    nb = get_info(contents, re_nocc_beta)

    if na != nb:
        sys.exit('numers of alpha and beta electron are not equal.\nexiting..')

    input_dict = {
            "title"         : "'"+os.path.splitext(logfile)[0]+"'",
            "iprint"         : '2',
            "gdictnfile"     : "'"+os.path.splitext(logfile)[0]+".F10'",
            "gintegfile"     : "'"+os.path.splitext(logfile)[0]+".F09'",
            "gbasisfile"     : "'"+os.path.splitext(logfile)[0]+".basinfo'",
            "nmos"           : (na+nb)/2,
            "maxiters"       : 100,
            "tstthr"         : '1.0d-10',
            "thresh"         : '1.0d-6',
            "alpha"          : '1.373',
            "beta"           : '2.0',
            "gamma"          : '2.0',
            "nppr"           : '3',
            "scfdamp"        : '0.5d0',
            "lrfun"          : '.false.',
            "lsym"           : '.false.',
            "lintsm"         : '.false.',
            }

    return input_dict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("logfile", help="gamess-us log file")
    parser.add_argument("-o", "--output", default="",
                        help="script output / dsfun input file")
    args = parser.parse_args()

    input_dct = get_input_data(args.logfile)

    if args.output == "":
        input_file = os.path.splitext(args.logfile)[0]+'_dsfun.inp'
    else:
        input_file = args.output

    write_input(input_file, input_dct)

if __name__ =="__main__":
    main()
