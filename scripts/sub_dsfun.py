#!/usr/bin/env python
import os
import re
import string
import subprocess
import sys
from optparse import OptionParser

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


def get_info(lines, rege):
    '''Get the infromation form contents based on the regular expression.'''

    cp = re.compile(rege)

    for line in lines:
        match = cp.search(line)
        if match:
            return match.group(1)

def main():
    p = OptionParser()
    p.add_option("-j",action="store",type="string",dest="jobname",default="")
    p.add_option("-t",action="store",type="string",dest="timeLimit",default="120:00:00")
    p.add_option("-m",action="store_true",dest="sendemail",default=False)
    p.add_option("-d",action="store_true",dest="norun",default=False)
    (opts, args) = p.parse_args()

    if opts.jobname == "":
        sys.exit("no jobname given. bye bye...")

    executable = '/home/lmentel/Source/vconstr/Source/bin/dsfun.x'
    workdir = os.getcwd()
    home    = os.getenv("HOME")
    scratch = os.getenv("TMPDIR")
    jobBase = os.path.splitext(opts.jobname)[0]
    output  = jobBase+'.out'
    error   = jobBase+'.err'
# get the gamess filenames from the dmft input
    f = open(opts.jobname, 'r')
    contents = f.readlines()
    f.close
    re_basis  = r"\s*gbasisfile=\'(.+)\'"
    re_dict   = r"\s*gdictnfile=\'(.+)\'"
    re_twoint = r"\s*gintegfile=\'(.+)\'"
    re_grid   = r"\s*gridfile=\'(.+)\'"
    basisfilename  = get_info(contents, re_basis)
    dictfilename   = get_info(contents, re_dict)
    twointfilename = get_info(contents, re_twoint)
    gridfilename   = get_info(contents, re_grid)
    points = 'points'

# open and  write the script PBS script file 
    scriptName = "run." + opts.jobname
    script = open(scriptName, 'w')

    script.write("#PBS -S /bin/bash\n")
    script.write("#PBS -l nodes=1:ppn=1\n")
    script.write("#PBS -l walltime="+opts.timeLimit+"\n\n")
#    script.write("#PBS -o {0:>s}.$PBS_O_JOBID \n".format(output))
#    script.write("#PBS -e {0:>s}.$PBS_O_JOBID \n\n".format(error))
    script.write("cd $PBS_O_WORKDIR\n")
    script.write("mkdir -p $TMPDIR/{0}\n".format(jobBase))
    script.write("cp {0} $TMPDIR/{1}/ \n".format(opts.jobname,  jobBase))
    script.write("cp {0} $TMPDIR/{1}/ \n".format(basisfilename, jobBase))
    script.write("cp {0} $TMPDIR/{1}/ \n".format(dictfilename,  jobBase))
    script.write("cp {0} $TMPDIR/{1}/ \n".format(twointfilename,jobBase))
    script.write("cp {0} $TMPDIR/{1}/ \n".format(gridfilename, jobBase))
    script.write("cd $TMPDIR/{0}\n\n".format(jobBase))  
    if opts.sendemail:
        script.write('''echo -e "Job STARTED\
        \t\\nPBS_JOBID:  $PBS_JOBID \
        \t\\nINPUT NAME: {0}\
        \t\\nWORKDIR:    {1} \
        \t\\nstarted at: `date`" | mail $USER -s "Job $PBS_JOBID"\n'''.format(opts.jobname,workdir))
    script.write("module load fortran/intel/64 c/intel/64 mkl/64 \n")
    script.write("{0:<s} {1:<s} \n".format(executable, opts.jobname))
    script.write("cp -t {0:<s} *.dat *.plt ".format(workdir))
    if opts.sendemail:
        script.write('''echo -e "Job FINISHED\
        \t\\nPBS_JOBID:  $PBS_JOBID \
        \t\\nINPUT NAME: {0}\
        \t\\nWORKDIR:    {1} \
        \t\\nfinished at: `date`" | mail $USER -s "Job $PBS_JOBID"'''.format(opts.jobname,workdir))
    script.close()

# submit the job to the queue if requested  
    if opts.norun:
        print "created job script: {0}\n NOT submitting to the queue\nbye...".format(scriptName) 
    else:
        print "created job script: {0}\n submitting to the queue".format(scriptName) 
        sublog = open(opts.jobname.replace(".inp",".sublog"),'w')
        proc = subprocess.Popen(["qsub",scriptName],stdout=sublog,stderr=sublog)
        sublog.close()

if __name__ == "__main__":
    main()
