#!/usr/bin/env python3

"""
:Author: Giuseppe Giovanni Nardone
:Contact: giuseppegiovanni.nardone [@] burlo.trieste.it
:Date: *2024.10.28

:Description:

Wrapper script to run the minION variant calling pipeline. It has two modes:
1 - Command line execution:
	- Will execute the pipeline given certain flags

2 - Interactive execution:
	- it will run the pipeline asking for arguments interactively

"""

# Wrapper script to launch the workflow/setup the snakemake environment and everything else that comes to mind!

import argparse
import pandas as pd 
import os
from datetime import datetime
import subprocess
import sys
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

import uuid

def check_env():
    envs = subprocess.check_output("""conda info --envs | cut -d " " -f1""", shell = True, encoding='utf8').split("\n")
    if 'snakemake8' in envs:
        print("########################################")
        print('Snakemake v8 environment available')
        return True
    else:
        print("########################################")
        print('Snakemake v8 environment not present. Installing...')
        return False

def check_samplesheet(samplesheet):
    ss = pd.read_csv(samplesheet, sep=" ", header=0, dtype='object')
    cols = ['SAMPLE_ID', "RAW_DATA", "ANALYSIS", "BARCODE"]
    print("########################################")
    print('Checking samplesheet')
    print("########################################")
    if set(cols).issubset(ss.columns) and len(ss.columns) == 4:
        print('Samplesheet is ready to be used')
        print("########################################")
    else:
        return sys.exit(f"""########################################\nSomething is wrong with the samplesheet file.\nPlease check if the separator is blank space and if contains the following columns:\n{" ".join(cols)}\n########################################""")

def run_cl(args):
    print("########################################")
    print("##### minION_v2 workflow #####") 
    print("########################################\n")

    #ss = args.samplesheet
    config = args.config
    logname = args.logname
    outdir = args.outdir
    dryrun = args.dryrun

    if not check_env():
        base_snakemake_dir=os.path.dirname((os.path.realpath(__file__)))
        env_file = os.path.join(base_snakemake_dir, "minION_workflow_v2.yaml")
        try:
            subprocess.check_output(f"""conda env create -f {env_file}""", shell = True, encoding = 'utf8', stderr=subprocess.STDOUT)
            print("Snakemake v8 installation completed")
        except subprocess.CalledProcessError as e:
            exit_code = e.returncode
            sys.exit(f"Something went wrong during the snakemake environment installation. Try to install it separately.\nExit code: {exit_code}")
    with open(config, 'r') as conf:
        c = load(conf, Loader=Loader)
    ss = c['paths']['samplesheet']
    check_samplesheet(samplesheet=ss)
    logdir = os.path.join(outdir, "Logs")
    if not os.path.isdir(logdir):
        os.makedirs(logdir)
    conda = subprocess.check_output("which conda", shell = True, encoding = 'utf8').strip("\n")
    base_snakemake_dir=os.path.dirname((os.path.realpath(__file__)))
    snakefile = os.path.join(base_snakemake_dir, "workflow", "Snakefile")
    profile = os.path.join(base_snakemake_dir, "workflow", "profile")
    activation_cmd = f"""eval "$({conda} shell.bash hook)" && {conda} activate snakemake8"""
    if not dryrun:
        snakemake_cmd = f"""snakemake -p -s {snakefile} --configfile {config} --directory {outdir} --workflow-profile {profile} --sdm conda env-modules --use-envmodules --use-conda --default-resources disk_mb=20000 1> {logdir}/{logname}.log 2> {logdir}/{logname}.err"""
    else:
        snakemake_cmd = f"""snakemake -n -p -s {snakefile} --configfile {config} --directory {outdir} --workflow-profile {profile} --sdm conda env-modules --use-envmodules --use-conda --default-resources disk_mb=20000 1> {logdir}/{logname}.log 2> {logdir}/{logname}.err"""
    #print(f"Executing...\n{snakemake_cmd}")
    #final_cmd = f"ml load singularity && {activation_cmd} && {snakemake_cmd}"
    #subprocess.call(final_cmd, shell=True, executable='/bin/bash')
    shebang="#!/usr/bin/env bash"
    filename = f"{uuid.uuid4().hex}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.sh"
    with open(os.path.join(outdir, filename), 'w') as script:
        script.write(f"{shebang}\n\nml load singularity\n{snakemake_cmd}")
    print(f"Executing...\n{snakemake_cmd}")
    final_cmd = f"""conda run -n snakemake8 bash {os.path.join(outdir, filename)} && rm {os.path.join(outdir, filename)}"""
    subprocess.call(final_cmd, shell=True, executable='/bin/bash')


def run(args):
    print("##### minION_v2 workflow #####\n")
    
    #ss = input("Enter the samplesheet path: ")
    config = input("Enter the config file path: ")
    logname = input("Enter the prefix for log files: ")
    outdir = input("Enter the output directory path: ")
    dryrun = input("Do you want to perform a dry-run? ")

    if not check_env():
        base_snakemake_dir=os.path.dirname((os.path.realpath(__file__)))
        env_file = os.path.join(base_snakemake_dir, "minION_workflow_v2.yaml")
        try:
            subprocess.check_output(f"""conda env create -f {env_file}""", shell = True, encoding = 'utf8', stderr=subprocess.STDOUT)
            print("Snakemake v8 installation completed")
        except subprocess.CalledProcessError as e:
            exit_code = e.returncode
            sys.exit(f"Something went wrong during the snakemake environment installation. Try to install it separately.\nExit code: {exit_code}")
    with open(config, 'r') as conf:
        c = load(conf, Loader=Loader)
    ss = c['paths']['samplesheet']

    check_samplesheet(samplesheet=ss)
    logdir = os.path.join(outdir, "Logs")
    if not os.path.isdir(logdir):
        os.makedirs(logdir)
    conda = subprocess.check_output("which conda", shell = True, encoding = 'utf8').strip("\n")
    base_snakemake_dir=os.path.dirname((os.path.realpath(__file__)))
    snakefile = os.path.join(base_snakemake_dir, "workflow", "Snakefile")
    profile = os.path.join(base_snakemake_dir, "workflow", "profile")
    activation_cmd = f"""eval "$({conda} shell.bash hook)" && {conda} activate snakemake8"""
    if not "y" in dryrun.lower():
        snakemake_cmd = f"""snakemake -p -s {snakefile} --configfile {config} --directory {outdir} --workflow-profile {profile} --sdm conda env-modules --use-envmodules --use-conda --default-resources disk_mb=20000 1> {logdir}/{logname}.log 2> {logdir}/{logname}.err"""
    else:
        snakemake_cmd = f"""snakemake -n -p -s {snakefile} --configfile {config} --directory {outdir} --workflow-profile {profile} --sdm conda env-modules --use-envmodules --use-conda --default-resources disk_mb=20000 1> {logdir}/{logname}.log 2> {logdir}/{logname}.err"""
    shebang="#!/usr/bin/env bash"
    filename = f"{uuid.uuid4().hex}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.sh"
    with open(os.path.join(outdir, filename), 'w') as script:
        script.write(f"{shebang}\n\nml load singularity\n{snakemake_cmd}")
    print(f"Executing...\n{snakemake_cmd}")
    #final_cmd = f"ml load singularity && {activation_cmd} && {snakemake_cmd}"
    final_cmd = f"""conda run -n snakemake8 bash {os.path.join(outdir, filename)} && rm {os.path.join(outdir, filename)}"""
    subprocess.call(final_cmd, shell=True, executable='/bin/bash')

# Set the parent parser
parent_parser = argparse.ArgumentParser(description='Wrapper script to run the minION workflow v2!')
# Add a child subparser
mode_parser = parent_parser.add_subparsers()
# Mode configuration - command line
cl_p = mode_parser.add_parser('cl',formatter_class=argparse.ArgumentDefaultsHelpFormatter,description='Launch the workflow using the command line')
cl_p.set_defaults(func=run_cl)
#cl_p.add_argument('-s', "--samplesheet", help="Samplesheet to use in the analyses", type=str, required=True)
cl_p.add_argument('-c', '--config', help='Config file to use', required=True, type=str)
cl_p.add_argument('-l', '--logname', help='Prefix of the log files', type=str, required=False, default=datetime.now().strftime('minION_v2_%Y-%m-%d_%H-%M-%S'))
cl_p.add_argument('-o', '--outdir', help='Name of the output directory in which store analyses results', required=True, type=str)
cl_p.add_argument('-n', '--dryrun', help='Perform a dry run to check if everything is okay', action="store_true", required=False)

interactive_p = mode_parser.add_parser('run',formatter_class=argparse.ArgumentDefaultsHelpFormatter,description='Launch the workflow interactively')
interactive_p.set_defaults(func=run)

# force the help message if the runner is called without arguments
if sys.argv[1] != 'run' and len(sys.argv) <= 2:
    sys.argv.append('--help')

options = parent_parser.parse_args()

options.func(options)
