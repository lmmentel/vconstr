#!/usr/bin/env python

import argparse
import os
import socket
import subprocess

def remove_dat(path, datfile):
    '''Remove the dat file from the ASCII scratch.'''
    if os.path.exists(os.path.join(path, datfile)):
        os.remove(os.path.join(path, datfile))

def main():
    '''
    Script for submitting gamessus jobs to the queue.
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument("input",
                        help='gamessus input file to be executed')
    parser.add_argument("-d",
                        "--dryrun",
                        action="store_true",
                        help="write the pbs script but don't submit to the queue, default-False")
    parser.add_argument("-t",
                        "--walltime",
                        default="120:00:00",
                        help="walltime in the format HH:MM:SS,\
                        default=120:00:00")
    args = vars(parser.parse_args())

def set_defaults(args):

    args['extrafiles'] = []
    args['gmsver'] = '00'
    args['nodes'] = '1'
    args['HOST'] = ''
    args['ppn'] = '1'
    args['usescratch'] = False
    args['queue'] = 'default'

    args['workdir'] = os.getcwd()
    args['home'] = os.getenv("HOME")
    if socket.gethostname() in ["login1.lisa.surfsara.nl", "login2.lisa.surfsara.nl"]:
        args['scratch'] = os.getenv("TMPDIR")
        args['rungms'] = os.getenv("RUNGMS")
        args['gmsver'] = os.getenv("GMSVER")
    args['local_scr'] = os.path.join(os.getenv('HOME'), 'scratch')
    args['jobname'] = os.path.splitext(args["input"])[0]
    args['outfile'] = args['jobname'] + ".log"
    args['errfile'] = args['jobname'] + ".err"
    args['datfile'] = args['jobname'] + ".dat"
    args['script_name'] = "run." + args['jobname']
    return args

def submit_pbs(args):
    '''
    Write the run script for PBS and submit it to the queue.
    '''

    args = set_defaults(args)
    remove_dat(args["local_scr"], args["datfile"])

    with open(args['script_name'], 'w') as script:
        script.write("#PBS -S /bin/bash\n")
        if args['HOST'] != "":
            script.write("#PBS -l nodes={0}:ppn={1}\n".format(args['HOST'], args['ppn']))
        else:
            script.write("#PBS -l nodes={0}:ppn={1}\n".format(args['nodes'], args['ppn']))
        if "mem" in args.keys():
            script.write("#PBS -l mem={0}\n".format(args["mem"]))
        script.write("#PBS -l walltime={0}\n\n".format(args['walltime']))
        #script.write('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/lib\n\n')
        script.write("#PBS -N {}\n".format(args["jobname"]))
        if args["queue"] != "default":
            script.write("#PBS -q {}\n".format(args["queue"]))
        script.write("cd $PBS_O_WORKDIR\n")
        if args['usescratch']:
            wrkdir = os.path.join(args['scratch'], args['jobname'])
            script.write("mkdir -p {}\n".format(wrkdir))
            files = args['jobname']
            if args['extrafiles']:
                files += ' ' + ' '.join(args['extrafiles'])
            script.write('cp -t {0} {1}\n'.format(wrkdir, files))
            script.write('cd {0}\n'.format(wrkdir))
        script.write("\n{0:<s} {1} {2:<s} {3:<s} \n".format(args['rungms'],
                     args['jobname'], args['gmsver'], args['ppn']))
        script.write("module load fortran/intel/64 c/intel/64 mkl/64 \n")
        script.write("\n{0:<s} {1:<s} {2:<s} {3:<s} \n".format(args['rungms'],
                     args['jobname'], args['gmsver'], args['ppn']))
        script.write("cp -t {0:<s} *.F09 *.F10 *.basinfo ".format(args['workdir']))

    # submit the job to the queue if requested
    if args['dryrun']:
        print("Created job script: {0}\n NOT submitting to the queue\nbye...".format(args['script_name']))
    else:
        print("Created job script: {0}\nsubmitting to the queue".format(args['script_name']))
        output = subprocess.check_output(["qsub", args['script_name']])
        pid = output.split(".")[0]
        return args['jobname'] + ".o" + pid

def submit_slurm(args):
    '''
    Write the run script for SLURM and submit it to the queue.
    '''

    with open(args['script_name'], 'w') as script:

        script.write("#!/bin/bash\n")
        script.write("#SBATCH -t {0:<s} \n".format(args["walltime"]))
        script.write("#SBATCH -N {0:<s} \n".format(args["nodes"]))
        script.write("#SBATCH --ntasks-per-node={0:<s} \n".format(args["ppn"]))
        script.write("#SBATCH -C {0:<s} \n".format(args["nodeType"]))
        script.write("#SBATCH -o {0:>s}.%J \n".format(args['outfile']))
        script.write("#SBATCH -e {0:>s}.%J \n".format(args['errfile']))

        if int(args["walltime"].split(':')[0])*60 + int(args["walltime"].split(':')[1]) > 60 and args["nodeType"] == "thin":
            script.write("#SBATCH -p normal \n\n")
        elif int(args["walltime"].split(':')[0])*60 + int(args["walltime"].split(':')[1]) < 60 and args["nodeType"] == "thin":
            script.write("#SBATCH -p short \n\n")
        elif args["nodeType"] == "fat":
            script.write("#SBATCH -p fat \n\n")

        if args["mail"] != "":
            script.write("#SBATCH --mail-type=ALL --mail-user={0:<s}\n\n".format(args["mail"]))
        script.write("mkdir -p $TMPDIR/{0}\n".format(args['jobname']))
        script.write("cp {0} $TMPDIR/{1} \n".format(args["input"], args['jobname']))
        script.write("cd $TMPDIR/{0}\n\n".format(args['jobname']))

        script.write("module load fortran/intel/64 c/intel/64 mkl/64 \n")
        script.write("{0:<s} {1:<s} 00 1 \n".format(args['rungms'], args['jobname']))
        script.write("cp -t {0:<s} *.F09 *.F10 *.basinfo ".format(args['workdir']))

    # submit the job to the queue if requested
    if args["dryrun"]:
        print("Created job script: {0}\n NOT submitting to the queue\nbye...".format(args['script_name']))
    else:
        print("Created job script: {0}\nsubmitting to the queue".format(args['script_name']))
        sublog = open(args['jobname'] + ".sublog", 'w')
        proc = subprocess.Popen(["sbatch", args['script_name']], stdout=sublog, stderr=sublog)
        sublog.close()

if __name__ == "__main__":
    main()
