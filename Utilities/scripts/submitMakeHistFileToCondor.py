#!/usr/bin/env python3

import argparse
import makeFileList 
from python import UserInput,OutputTools
from python import ConfigureJobs
import datetime
import logging
import os
import shutil
import glob
import tarfile
import math
import re
import subprocess
import logging

def getComLineArgs():
    parser = UserInput.getDefaultParser(False)
    parser.add_argument("-d", "--submit_dir", type=str,
                        required=True, help="Output directory")
    parser.add_argument("--debug", action='store_true',
        help="Print verbose info")
    parser.add_argument("-n", "--files_per_job", type=int,
                        default=3, help="Number of files per job")
    parser.add_argument("--submit", action='store_true',
                        help="Just make the directory, don't submit")
    parser.add_argument("--local", action='store_true',
                        help="Use local files, e.g., file_path")
    parser.add_argument("--input_tier", type=str,
        default="", help="Selection stage of input files")
    parser.add_argument("--memory", type=int,
        default=2000, help="Memory to request for condor jobs")
    parser.add_argument("-q", "--queue", type=str,
        default="longlunch", help="lxplus queue, or 'uw' for wisconsin settings")
    parser.add_argument("-m", "--merge", nargs=2, type=str,
        metavar=("mergedFileName", "completeFraction"),
        default=None, help="Merge outputs from all jobs to file (submit as DAG)" \
        "if > completeFraction (in %%) jobs complete")
    parser.add_argument("--numFiles", default=-1, type=int, help="Number of files to process")
    parser.add_argument("-j", "--numCores", default=1, type=int, help="Number of cores for each job")
    parser.add_argument("--force", action='store_true',
        help="Force overwrite of existing directories")
    parser.add_argument("--removeUnmerged", action='store_true',
        help="Remove unmerged ROOT files (requires DAG)")

    args = parser.parse_args()
    if args.removeUnmerged and not args.merge:
        parser.error("--removeUnmerged requires --merge")

    return vars(args)

def makeSubmitDir(submit_dir, force):
    log_dir = submit_dir + "/logs"
    if os.path.isdir(submit_dir) and force:
        logging.warning("Overwriting directory %s" % submit_dir)
        shutil.rmtree(submit_dir)
    elif os.path.isdir(submit_dir):
       raise IOError("Submit directory %s already exists! Use --force to overrite." % submit_dir) 
    os.makedirs(log_dir)

def setupMergeStep(submit_dir, queue, numjobs, merge, removeUnmerged):

    merge_file = merge[0]
    completeFraction = float(merge[1])
    if completeFraction < 0 or completeFraction > 1.:
        raise InvalidArgument("completeFraction must be between 0 and 1, found %f" % completeFraction)

    template_dict = {
        "queue" : queue,
        "merge_file" : merge_file
    }
    template = "Templates/CondorSubmit/merge_template.jdl"
    outfile = "/".join([submit_dir, "merge.jdl"])
    ConfigureJobs.fillTemplatedFile([template], outfile, template_dict)

    template = "Templates/CondorSubmit/merge.sh"
    outfile = "/".join([submit_dir, "merge.sh"])
    ConfigureJobs.fillTemplatedFile([template], outfile, {"CMSSW_RELEASE_BASE" : os.environ["CMSSW_BASE"]})

    template = "Templates/CondorSubmit/submit_and_merge_template.dag"
    outfile = "/".join([submit_dir, "submit_and_merge.dag"])
    ConfigureJobs.fillTemplatedFile([template], outfile, 
            {"minComplete" : int(completeFraction*numjobs),
            "postMerge" : ("SCRIPT POST B removeRootFiles.sh %s" % merge_file) if removeUnmerged else ""})

    for f in ["list_infiles.sh", "completed.sh", "removeRootFiles.sh"]:
        shutil.copy("Templates/CondorSubmit/%s" % f, "/".join([submit_dir, f]))

def copyLibs():
    libdir = "%s/src/lib" % os.environ["CMSSW_BASE"]
    if os.path.isdir(libdir):
        shutil.rmtree(libdir)
    os.mkdir(libdir)
    
    cmssw_libdir = "/".join([os.environ["CMSSW_BASE"], "lib", os.environ["SCRAM_ARCH"], "*SelectorTools*"])
    for i in glob.glob(cmssw_libdir):
        shutil.copyfile(i, '/'.join([libdir, os.path.basename(i)]))

# Needed on lxplus and uwlogin, where the afs permissions are set
# very tight and don't let condor access some dumb file it needs
# TODO: Understand what it needs and why
def modifyAFSPermissions():
    if not re.findall("system:anyuser *.*l", subprocess.check_output(["fs", "la"], encoding='UTF-8')):
        subprocess.call(["find", os.environ["CMSSW_BASE"], "-type", "d",
            "-exec", "fs", "setacl", "-dir", "{}", "-acl", "system:anyuser", "rl", ";"])
        raise OSError("AFS permissions have been relaxed for condor submission. You should recompile and resubmit")

def copyDatasetManagerFiles(analysis):
    manager_name = "AnalysisDatasetManager"
    manager_path = ConfigureJobs.getManagerPath()

    cwd = os.getcwd()
    os.chdir("%s/src" % os.environ["CMSSW_BASE"])
    if os.path.isdir(manager_name):
        shutil.rmtree(manager_name)

    info_dir = manager_name+"/FileInfo"
    plot_dir = manager_name+"/PlotObjects"
    os.makedirs(info_dir)
    os.makedirs(plot_dir)

    for paths in [[manager_path, manager_name, "Utilities"],
                    [manager_path, info_dir, "data"],
                    [manager_path, info_dir, "montecarlo"],
                    [manager_path, plot_dir, analysis]]:
        path = '/'.join(paths)
        if os.path.isdir(path):
            shutil.copytree(path, '/'.join(paths[1:]))
        else:
            for d in glob.glob(path+'*'):
                shutil.copy(d, d.replace(manager_path+"/", ""))

    os.chdir(cwd)

# TODO: Check if this is needed at UW. I think it isn't
def copyGridCertificate():
    proxypath = "/tmp/x509up_u%s" % os.getuid() if not \
                ("X509_USER_PROXY" in os.environ and os.path.isfile(os.environ["X509_USER_PROXY"])) \
            else os.environ["X509_USER_PROXY"] 
    shutil.copy(proxypath, "userproxy")

def tarAnalysisInfo(condor_dir, tarball_name):
    tarname = condor_dir+"/"+tarball_name
    currd = os.getcwd()
    os.chdir("../../")
    base = "Analysis/SelectorTools"
    with tarfile.open(tarname, "w:gz") as tar:
        tar.add(f"lib")
        tar.add(f"AnalysisDatasetManager")
        tar.add(f"{base}/data")
        tar.add(f"{base}/Utilities")
        tar.add(f"{base}/interface")
        tar.add(f"{base}/src")
    shutil.rmtree("lib")
    shutil.rmtree("AnalysisDatasetManager")
    os.chdir(currd)

def getUWCondorSettings():
    return """# (wisconsin-specific) tell glideins to run job with access to cvmfs (via parrot)
        +RequiresCVMFS       = True
        +RequiresSharedFS    = True
        +IsFastQueueJob      = True
        Requirements         = TARGET.Arch == "X86_64" && IsSlowSlot=!=true && (MY.RequiresSharedFS=!=true || TARGET.HasAFS_OSG) && (TARGET.HasParrotCVMFS=?=true || (TARGET.UWCMS_CVMFS_Exists  && TARGET.CMS_CVMFS_Exists))
    """

def writeSubmitFile(submit_dir, analysis, selection, input_tier, queue, memory, filelist, numfiles, numCores, nPerJob, selArgs):
    if nPerJob < numCores:
        logging.warning("Number of cores is less than number of files per job. Setting instead to %i" % nPerJob)
        nPerJob

    extraArgs = "--debug --compress %s" % ("" if not selArgs else ("--selectorArgs "+" ".join(selArgs)))
    if numCores != 1:
        extraArgs = "-j %i %s" % (numCores, extraArgs)

    template_dict = {
        "analysis" : analysis,
        "selection" : selection,
        "input_tier" : input_tier,
        "queue" : queue,
        "memory" : memory,
        "filelist" : filelist.split(".txt")[0],
        "nPerJob" : nPerJob,
        "numCores" : numCores,
        "nJobs" : int(math.ceil(float(numfiles)/nPerJob)),
        "extraArgs" : extraArgs
    }

    template = "Templates/CondorSubmit/submit_template.jdl"
    outfile = "/".join([submit_dir, "submit.jdl"])
    ConfigureJobs.fillTemplatedFile([template], outfile, template_dict)

def writeWrapperFile(submit_dir, tarball_name):
    template_dict = { 
            "CMSSW_RELEASE_BASE" : os.environ["CMSSW_RELEASE_BASE"],
            "tarball" : tarball_name
    }
    template = "Templates/CondorSubmit/wrapRunSelector.sh"
    outfile = "/".join([submit_dir, "wrapRunSelector.sh"])
    ConfigureJobs.fillTemplatedFile([template], outfile, template_dict)

def writeMetaInfo(submit_dir, filename):
    with open("/".join([submit_dir, filename]), "w") as metafile:
        metafile.write(OutputTools.getScriptCall()+"\n\n")
        metafile.write("Script called at %s \n" % datetime.datetime.now())
        metafile.write("git hash: " + OutputTools.gitHash()+"\n")
        metafile.write("git diff: " + OutputTools.gitDiff()+"\n")

def submitDASFilesToCondor(filenames, submit_dir, analysis, selection, input_tier, queue, memory,
        numPerJob, force, das, selArgs, merge, removeUnmerged, maxFiles, numCores):
    makeSubmitDir(submit_dir, force)
    copyLibs()
    copyDatasetManagerFiles(analysis)
    if "/afs" in os.getcwd()[:4]:
        modifyAFSPermissions()

    filelist_name = '_'.join(filenames[:max(len(filenames), 4)])
    filelist_name = filelist_name.replace("*", "ALL")
    filelist = '/'.join([submit_dir, filelist_name+'_filelist.txt'])
    numfiles = makeFileList.makeFileList(filenames, filelist, analysis, input_tier, das)
    if maxFiles > 0 and maxFiles < numfiles:
        numfiles = maxFiles

    #TODO: I don't think there's any harm in addition the accounting group, but
    # it doesn't do anything if you aren't a member of CMST3 group
    cernk = "kelong" in os.getlogin()
    if queue == 'uw':
        getUWCondorSettings()
    elif queue == 'mit':
        #queue = 'requirements = HAS_CVMFS_cms_cern_ch'
        #queue = '+SingularityImage = "/cvmfs/cernvm-prod.cern.ch/cvm3"'
        queue = ('Requirements = regexp("T3BTCH*", MACHINE) \n'
                 '+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el7:latest"'
                )
        #queue = '+REQUIRED_OS = "rhel7"'
    else:
        queue = '+JobFlavour = "{0}"\n+AccountingGroup = "group_u_CMST3.all"'.format(queue) \

    writeSubmitFile(submit_dir, analysis, selection, input_tier, queue, memory, filelist_name, numfiles, numCores, numPerJob, selArgs)
    if merge:
        setupMergeStep(submit_dir, queue, math.ceil(numfiles/numPerJob), merge, removeUnmerged)

    tarball_name = '_'.join([analysis, "AnalysisCode.tgz"])
    writeWrapperFile(submit_dir, tarball_name)
    tarAnalysisInfo(submit_dir, tarball_name)
    writeMetaInfo(submit_dir, "metaInfo.txt")

def main():
    args = getComLineArgs()
    logging.basicConfig(level=(logging.DEBUG if args['debug'] else logging.INFO))
    submitDASFilesToCondor(args['filenames'], args['submit_dir'], args['analysis'], 
        args['selection'], args['input_tier'], args['queue'], args['memory'], args['files_per_job'], args['force'], 
        not args['local'], args['selectorArgs'], args['merge'], args['removeUnmerged'],
        args['numFiles'], args['numCores'])
    if args['submit']:
        command = 'condor_submit' if not args['merge'] else 'condor_submit_dag'
        submitfile = 'submit.jdl' if not args['merge'] else 'submit_and_merge.dag'
        os.chdir(args['submit_dir'])
        subprocess.call([command, '/'.join([args['submit_dir'], submitfile])])

if __name__ == "__main__":
    main()
