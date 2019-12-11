import os
import shutil
import sys
import time
import argparse
import re
import logging
from logging.handlers import RotatingFileHandler
from subprocess import call, Popen, PIPE, DEVNULL
from report import getReport
from datetime import datetime


'''
FLORA stands for (pipeline) For LOng Read Assembly
(subject to change)
'''

# Option parser
PARSER = argparse.ArgumentParser()
PARSER.add_argument("-i", "--contig", dest="contig", default=None, help="Input contig files.")
PARSER.add_argument("-r", "--read", dest="read", default=None, help="Input long-read in fastq format. WARNING, you must have the read in fasta format with the exact same name in the same directory than the fastq file.")
PARSER.add_argument("-c", "--correction", dest="correct", default=None, help="Input read used for correction")
PARSER.add_argument("-sr1", "--shortread1", dest="shortread1", default=None, help="Input of paired read used for hybrid assembly")
PARSER.add_argument("-sr2", "--shortread2", dest="shortread2", default=None, help="Input of paired read used for hybrid assembly")
PARSER.add_argument("-a", "--assembler", dest="assembler", default=None, choices=["Flye","F","WTDBG2","W","Spades","S","MaSuRCA","M","Canu","C"], help="ID of the assembler you want to use (option --assembler_list or -al for the list)")
PARSER.add_argument("-po", "--polisher", dest="polisher", default="Pilon", choices=["Pilon","P","Racon","R"], help="ID of the polisher you want to use (option --polisher_list or -pl for the list, default is Pilon)")
PARSER.add_argument("-p", "--pattern", dest="pattern", default=None, help="Launch tools in the order of the pattern (ex: FTAQPsQPsQ). See --tuto (-u) for a little help and list of letter available")
PARSER.add_argument("-u", "--tuto", dest="tuto", default=None, action="store_true", help="Display a little help and list of letter available for the pattern option")
PARSER.add_argument("-t", "--thread", dest="thread", default=8, type=int, help="Max number of CPU thread available for the tools (default = 8).")
PARSER.add_argument("-m", "--ram", dest="ram", default="8", help="Max giga of RAM available for the tools (default = 8 G).")
PARSER.add_argument("-q", "--quality", dest="quality", default=None, help="Quality threshold for filter by NanoFilt (only usefull if Filter step (F) is in pattern).")
PARSER.add_argument("-l", "--length", dest="length", default=None, help="Length threshold for filter by NanoFilt (only usefull if Filter step (F) is in pattern).")
PARSER.add_argument("-e", "--estimate", dest="estimate", default=None, help="Estimate size of the input genome (ex: 80m).")
PARSER.add_argument("-k", "--kmer", dest="kmer", default=25, type=int, help="K length of mer (default = 25).")
PARSER.add_argument("-al", "--assembler_list", dest="assembler_list", default=None, action="store_true", help="Display the list of ID of assembler available")
PARSER.add_argument("-pl", "--polisher_list", dest="polisher_list", default=None, action="store_true", help="Display the list of ID of polisher available")
PARSER.add_argument("-o", "--output", dest="output", default="run", help="ID of the run. It will use this ID for output prefix.")
PARSER.add_argument("-d", "--dir", dest="dir", default="~/FLORA_OUT/", help="Directory to stock result (default = ~/FLORA_OUT/).")
PARSER.add_argument("-re","--retry", dest="retry", default=None, help="Allow to retry a run where FLORA failed. Need a log file as input.")
PARSER.add_argument("-ro","--readopt", dest="readopt", default="nano-raw", choices=["nano-raw","nano-corr"], help="Read mode for Flye and Canu.")
PARSER.add_argument("-ao","--aligneropt", dest="aligneropt", default="minimap2", choices=["minimap2","bwa"], help="Aligner choice.")
PARSER.add_argument("-f", "--fullpath", dest="fullpath", default=False, action="store_true", help="Show full name of file in the html report. Usefull when using result from two or more run.")
store_action_dict=vars(PARSER)['_option_string_actions']

if __name__ == "__main__":

    #
    if len(sys.argv)==1 :
        print()
        PARSER.print_help(sys.stderr)
        print()
        sys.exit(1)
    args = PARSER.parse_args()
    #
    opt = list(store_action_dict.keys())
    dictopt = dict(zip(opt[::2], opt[1::2]))

# Formatter for log file 
    class MyFormatterFile(logging.Formatter):

        err_fmt = "\n%(asctime)s :: %(levelname)s :: %(message)s"
        dbg_fmt  = "\n%(asctime)s :: %(levelname)s :: %(message)s"
        info_fmt = "\n\n%(asctime)s :: %(levelname)s :: %(message)s"

        def __init__(self):
            super().__init__(fmt="%(levelno)d: %(msg)s", datefmt=None, style='%')  

        def format(self, record):
            format_orig = self._style._fmt
            if record.levelno == logging.DEBUG:
                self._style._fmt = MyFormatterFile.dbg_fmt
            elif record.levelno == logging.INFO:
                self._style._fmt = MyFormatterFile.info_fmt
            elif record.levelno == logging.ERROR:
                self._style._fmt = MyFormatterFile.err_fmt
            result = logging.Formatter.format(self, record)
            self._style._fmt = format_orig
            return result
    class MyFormatterStream(logging.Formatter):

        err_fmt = "\t---- %(levelname)s ----\n\n%(message)s\n"
        dbg_fmt  = "%(levelname)s :: %(message)s"
        info_fmt = "\t---- %(message)s ----\n" 

        def __init__(self):
            super().__init__(fmt="%(levelno)d: %(msg)s", datefmt=None, style='%')  

        def format(self, record):
            format_orig = self._style._fmt
            if record.levelno == logging.DEBUG:
                self._style._fmt = MyFormatterStream.dbg_fmt
            elif record.levelno == logging.INFO:
                self._style._fmt = MyFormatterStream.info_fmt
            elif record.levelno == logging.ERROR:
                self._style._fmt = MyFormatterStream.err_fmt
            result = logging.Formatter.format(self, record)
            self._style._fmt = format_orig
            return result
# Checking log file input if args.retry is set
    if args.retry :
        logFile = args.retry
        if not os.path.isfile(logFile) :
            sys.exit("ERROR, input of retry option doesn't exist, need a log file.")
        elif logFile.split(".")[-1] != "log" :
            sys.exit("ERROR, input of retry option isn't a log file.")
# Getting old options from log file 
        with open(logFile,"r") as log :
            for i in log :
                if "DEBUG" in i :
                    cmd=i.split(" :: ")[-1]
                    break
            res=re.findall("(-\w+|--\w+) ?((/?\w+)*(\.\w+)?)",cmd)
            if res: 
                for i in res :
                    ar=i[0].strip()
                    if ar[0:2] == "--" : 
                        if ar in dictopt.values():
                            op=ar[2:]
                            if i[1]=="":
                                val=True
                            else :
                                val=i[1]
                    else:
                        if ar in dictopt.keys():
                            op=dictopt[ar][2:]
                            if i[1]=="":
                                val=True
                            else :
                                val=i[1]
                    args.__setattr__(op,val)
    else :
        logFile = ""

    if args.dir[-1]!="/":
        args.dir+="/"
#    if args.dir[0]!="/":
#        if args.dir[0]!="~":
#            args.dir="~/"+args.dir
#        elif args.dir[1]!="/":
#            args.dir="~/"+args.dir[1:]

    if args.ram[-1]=="G" :
        args.ram=args.ram[:-1]

# Setting up the variable step which is used for the log file
    global step 
    step=1
# List of assembler available
    assemblerList = ["Flye","F","Canu","C","WTDBG2","W","Masurca","M"]
    if args.assembler_list: 
        i = 0
        print("\n\tASSEMBLER LIST\n")
        while i < len(assemblerList)-1:
            print('{:.<30s}{:>1s}'.format(assemblerList[i],assemblerList[i+1]))
            i += 2
        print()
        sys.exit()

# List of polisher available
    polisherList = ["Pilon","P","Racon","R"]
    if args.polisher_list: 
        i = 0
        print("\n\tPOLISHER LIST\n")
        while i < len(polisherList)-1:
            print('{:.<30s}{:>1s}'.format(polisherList[i],polisherList[i+1]))
            i += 2
        print()
        sys.exit()

# List of letter accepted for the pattern option and their meaning (used for the tutorial)
    letterList = ["E","R","T","F","C","A","Q","Ps","Pl"]
    patternList = ["Estimate","Read quality","Trim","Filter","Correct","Assemble","Quality control","Polish short-read","Polish long-read"]

    if args.tuto:
        i = 0
        print("\n\tTUTORIAL")
        print("\nFLORA need a pattern to known which tools to launch when.\nFor example the pattern AQ launch a Assembly then Quality (Quast and Busco).\n\nYou can chain multiple times the same tool.\nFor example, the pattern FAPlQPsQ launch these step in this order : \n\n\t- Filter \n\t- Assembly \n\t- Polish with long-read \n\t- Quality \n\t- Polish with short-read \n\t- Quality \n\nSee list below for available option :")
        print()
        while i < len(patternList):
            print('{:.<30s}{:>1s}'.format(patternList[i],letterList[i]))
            i += 1
        print()
        sys.exit()

    prefix = args.output+"_FLORA"

# Check if specific reads for correction are given, if not, FLORA use the long-read for self-correction / polishing
    if not args.correct :
        args.correct = args.read

# Control of read format which is stored in variable typ  
    if args.read:
        ext = args.read.split(".")[-1]
        if ext in ["fa", "fasta","fas"]:
            typ = "fasta"
        elif ext in ["fq", "fastq"]:
            typ = "fastq"
        else:
            print(args.read)
            sys.exit("ERROR, wrong format for read input (fasta or fastq only).")

# Control of pattern validity 
    if args.pattern:
        if args.pattern[0].islower():
            sys.exit("ERROR, lowercase char at beginning of pattern option.")
        i = 0
        pat = []
        while i < len(args.pattern):
            if args.pattern[i].islower():
                pat[-1] = pat[-1]+args.pattern[i]
            else:
                pat.append(args.pattern[i])
            i += 1
        wrongChar = ""
        for i in pat:
            if i not in letterList:
                if wrongChar:
                    wrongChar = wrongChar+" - "+i
                else:
                    wrongChar = wrongChar+i
        if wrongChar:
            sys.exit("ERROR, wrong char in pattern option ("+wrongChar+").")
    else:
        sys.exit("ERROR, need a pattern (see --tuto).")

# Comment the following 4 lines if you want to do multiple filter and/or trim step
    if pat.count("T") > 1:
        sys.exit("ERROR, can't do more than one trim step (T).")
    if pat.count("F") > 1:
        sys.exit("ERROR, can't do more than one filter step (F).")

    cwd = os.path.dirname(os.path.abspath(__file__))
# Checking if config file exists
    if not os.path.isfile(cwd+"/config") :
        sys.exit("ERROR, missing config file, generate it with python3 setup.py.")

# Config file store all software path
    with open(cwd+"/config","r") as conf :
        for i in conf :
            if "flye" in i:
                Flye_path = i.split(" = ")[-1][:-1]
            elif "porechop" in i:
                Porechop_path = i.split(" = ")[-1][:-1]
            elif "wtdbg2" in i:
                Wtdbg2_path = i.split(" = ")[-1][:-1]
            elif "NanoFilt" in i:
                NanoFilt_path = i.split(" = ")[-1][:-1]
            elif "NanoPlot" in i:
                NanoPlot_path = i.split(" = ")[-1][:-1]
            elif "ropebwt2" in i:
                Ropebwt2_path = i.split(" = ")[-1][:-1]
            elif "quast" in i:
                Quast_path = i.split(" = ")[-1][:-1]
            elif "busco" in i:
                Busco_path = i.split(" = ")[-1][:-1]
            elif "racon" in i:
                Racon_path = i.split(" = ")[-1][:-1]
            elif "pilon" in i:
                Pilon_path = i.split(" = ")[-1][:-1]
            elif "fmlrc" in i:
                Fmlrc_path = i.split(" = ")[-1][:-1]
            elif "minimap2" in i:
                Minimap2_path = i.split(" = ")[-1][:-1]
            elif "samtools" in i:
                Samtools_path = i.split(" = ")[-1][:-1]
            elif "masurca" in i:
                Masurca_path = i.split(" = ")[-1][:-1]
            elif "canu" in i:
                Canu_path = i.split(" = ")[-1][:-1]
            elif "bwa" in i:
                Bwa_path = i.split(" = ")[-1][:-1]
            elif "jellyfish" in i:
                Jellyfish_path = i.split(" = ")[-1][:-1]

# Threads split for some parallel task
    gt=round(int(args.thread)/4*3)
    pt=round(int(args.thread)/4)
# Timestamp setup
    dateTimeObj = datetime.now()
    timestamp = dateTimeObj.strftime("%d-%m-%y_%Hh%Mm%Ss")

def flye(read,did):
    global step
    pdir=args.dir+str(did)+"-flye/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    logger.info(str(step)+". STARTING ASSEMBLY USING FLYE")
    time.sleep(2)
    step+=1
    contigR = pdir+"assembly.fasta"
    if args.retry :
        if not os.path.isfile(contigR):
            fly=Popen([Flye_path+"/flye","--"+args.readopt,read,"-t",str(args.thread),"--genome-size",args.estimate,"--out-dir",pdir])
            fly.wait()
            cmd=" ".join(fly.args)
            logger.debug(cmd)
    else :
        fly=Popen([Flye_path+"/flye","--"+args.readopt,read,"-t",str(args.thread),"--genome-size",args.estimate,"--out-dir",pdir])
        fly.wait()
        cmd=" ".join(fly.args)
        logger.debug(cmd)
#    if fly.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
    Ndid = did + 1
    return (contigR,Ndid)

def canu(read,did):
    global step
    pdir=args.dir+str(did)+"-canu/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    logger.info(str(step)+". STARTING ASSEMBLY USING CANU")
    time.sleep(2)
    step+=1
    if args.readopt=="nano-raw":
        canuopt="nanopore-raw"
    else :
        canuopt="nanopore-corr"
    contigR = pdir+prefix+".contigs.fasta"
    if args.retry :
        if not os.path.isfile(contigR):
            can=Popen([Canu_path+"/canu","-"+canuopt,read,"genomeSize="+args.estimate,"-d",pdir,"-p",prefix,"useGrid=false"])
            can.wait()
            cmd=" ".join(can.args)
            logger.debug(cmd)
    else :
        can=Popen([Canu_path+"/canu","-"+canuopt,read,"genomeSize="+args.estimate,"-d",pdir,"-p",prefix,"useGrid=false"])
        can.wait()
        cmd=" ".join(can.args)
        logger.debug(cmd)
#    if fly.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
    Ndid = did + 1
    return (contigR,Ndid)


def wtdbg2(read,did):
    global step
    pdir=args.dir+str(did)+"-wtdbg2/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    outP=pdir+prefix
    logger.info(str(step)+". STARTING ASSEMBLY USING WTDBG2")
    time.sleep(2)
    step+=1
    if args.retry :
        if not os.path.isfile(outP+".ctg.fa"):
            wtd=Popen([Wtdbg2_path+"/wtdbg2","-t"+str(args.thread),"-x","ont","-g",args.estimate,"-i",read,"-fo",outP])
            wtd.wait()
            cmd=" ".join(wtd.args)
            logger.debug(cmd) 
#    if wtd.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
            cns=Popen([Wtdbg2_path+"/wtpoa-cns","-t",str(args.thread),"-i",outP+".ctg.lay.gz","-fo",outP+".ctg.fa"])
            cns.wait()
            cmd=" ".join(cns.args)
            logger.debug(cmd)
    else :
        wtd=Popen([Wtdbg2_path+"/wtdbg2","-t"+str(args.thread),"-x","ont","-g",args.estimate,"-i",read,"-fo",outP])
        wtd.wait()
        cmd=" ".join(wtd.args)
        logger.debug(cmd) 
#    if wtd.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
        cns=Popen([Wtdbg2_path+"/wtpoa-cns","-t",str(args.thread),"-i",outP+".ctg.lay.gz","-fo",outP+".ctg.fa"])
        cns.wait()
        cmd=" ".join(cns.args)
        logger.debug(cmd)

#    if cns.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
    Ndid = did + 1 
    return(outP+".ctg.fa",Ndid)


#def spades(sRead1,lRead,sRead2=""):
#    print(sRead1)
#    print(sRead2)
#    print(lRead)
#    return

def masurca(sRead1,sRead2,lRead,did):
    write=False
    global step
    pdir=args.dir+str(did)+"-masurca/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    logger.info(str(step)+". STARTING HYBRID ASSEMBLY USING MASURCA")
    time.sleep(2)
    step+=1
    config_path=pdir+"example_config"
    mas=Popen([Masurca_path+"/masurca","-g",pdir+"example_config"])
    mas.wait()
    if not os.path.isfile(config_path):
        sys.exit("ERROR, no example of config file found.")
    with open(config_path,"r") as exconf, open(pdir+"mas_config.txt","w") as masconf :
        for i in exconf :
            if i[0]=="#" or i[0]=="\n":
                continue
            if i[0:4] in ["JUMP","PACB","OTHE"]:
                continue
            if i[0:2]=="PE" :
                if sRead2!="" :
                    sRead2=os.path.abspath(sRead2)
                masconf.write("PE= sr 150 20 "+os.path.abspath(sRead1)+" "+sRead2+"\n")
                masconf.write("NANOPORE="+os.path.abspath(lRead)+"\n")
                continue
            if i[0:5]=="NUM_T" :
                masconf.write("NUM_THREADS = "+str(args.thread)+"\n")
                continue
            masconf.write(i)
    if os.path.isfile(config_path):
        os.remove(config_path)
    mas=Popen([Masurca_path+"/masurca",pdir+"mas_config.txt","-o",pdir+"assemble.sh"])
    mas.wait()
    cmd=" ".join(mas.args)
    logger.debug(cmd)
    print()
    tcwd=os.getcwd()
    os.chdir(pdir)
    asse=Popen(["./assemble.sh"])
    asse.wait()
    cmd=" ".join(asse.args)
    logger.debug(cmd)
    os.chdir(tcwd)
    Ndid = did + 1 
    for i in os.listdir(pdir) :
        if os.path.isdir(i) and i[0:2]=="CA":
            contig_out = pdir+i+"/final.genome.scf.fasta"
    print(contig_out) 
    return(contig_out,Ndid)

def racon(mode,times,contig,read,did):
    global step
    pdir=args.dir+str(did)+"-racon/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    logger.info(str(step)+". STARTING ALIGNMENT OF CONTIG FILE AND READS")
    time.sleep(2)
    step+=1
    samout = pdir+prefix+"_minimap2_"+times+".sam"
    with open(samout,"w") as alout:
        if mode=="correction":
            if args.aligneropt == "bwa" :
                ind=Popen([Bwa_path+"/bwa","index",contig])
                ind.wait()
                mini=Popen([Bwa_path+"/bwa","mem","-t", str(args.thread), contig, read], stdout=PIPE)
            else :
                mini=Popen([Minimap2_path+"/minimap2","-x","sr","-t",str(gt),"-L","-a",contig, read], stdout=alout,stderr=PIPE)
        else:
            mini=Popen([Minimap2_path+"/minimap2","-t",str(gt),"-L","-a",contig, read], stdout=alout,stderr=PIPE)
        err=mini.communicate()
        cmd=" ".join(mini.args)+" > "+samout
        logger.debug(cmd)
        if mini.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
    logger.info(str(step)+". STARTING RACON "+mode.upper())
    step+=1
    time.sleep(2)
#    if mode == "polish" :
#        patte=re.compile(".*(_polished_x(\d+)).*")
#        mat=patte.match(nameC)
#        if mat :
#            nameC=nameC.replace(mat.group(1),"_polished_x"+str(int(mat.group(2))+1))
#    elif mode == "correction":
#        patte=re.compile(".*(_corrected_x(\d+)).*")
#        mat=patte.match(nameC)
#        if mat :
#            nameC=nameC.replace(mat.group(1),"_corrected_x"+str(int(mat.group(2))+1))
    con_r=pdir+args.output+"_racon"+times+".fasta"
    with open(con_r,"w") as conout :
        rac = Popen([Racon_path+"/racon", "-t",str(args.thread),read,samout,contig], stdout=conout,stderr=PIPE)
    err=rac.communicate()
    print()
    cmd=" ".join(rac.args)+" > "+con_r
    logger.debug(cmd)
    if rac.returncode != 0 :
        logger.error(err[-1].decode("utf-8"))
        sys.exit(1)
    Ndid = did + 1
    return (con_r,Ndid)


def pilon(mode,times,contig,read,did):
    global step
    pdir=args.dir+str(did)+"-pilon/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    logger.info(str(step)+". STARTING ALIGNMENT OF CONTIG FILE AND READS")
    time.sleep(2)
    step+=1
    bamout = pdir+prefix+"_minimap2_"+times+".bam"
    sortedout = pdir+prefix+"_minimap2_sorted_"+times+".bam"
    if not os.path.isfile(sortedout):
        with open(bamout,"w") as alout:
            if mode=="correction":
                if args.aligneropt == "bwa" :
                    ind=Popen([Bwa_path+"/bwa","index",contig])
                    ind.wait()
                    mini=Popen(Bwa_path+["/bwa","mem","-t", str(args.thread), contig, read], stdout=PIPE)
                else :
                    mini=Popen([Minimap2_path+"/minimap2","-x","sr","-t",str(gt),"-L","-a",contig, read], stdout=PIPE)
            else:
                mini=Popen([Minimap2_path+"/minimap2","-t",str(gt),"-L","-a",contig, read], stdout=PIPE)
            samt=Popen([Samtools_path+"/samtools","view","-@"+str(pt),"-Sb","-"],stdin=mini.stdout,stdout=alout,stderr=PIPE)
            err=samt.communicate()
            cmd=" ".join(mini.args)+" | "+" ".join(samt.args)
            logger.debug(cmd)
            if samt.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
            sams=Popen([Samtools_path+"/samtools","sort","-@"+str(args.thread),"-o",sortedout,bamout],stderr=PIPE)
            err=sams.communicate()
            cmd=" ".join(sams.args)
            logger.debug(cmd)
            if sams.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
            ind=Popen([Samtools_path+"/samtools", "index","-@"+str(args.thread), sortedout],stderr=PIPE)
            err=ind.communicate()
            cmd=" ".join(ind.args)
            logger.debug(cmd)
            if ind.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
    Ncontig="outPilon_"+times
    Ncontig_return=pdir+Ncontig+".fasta"
    logger.info(str(step)+". STARTING PILON "+mode.upper())
    step+=1
    time.sleep(2)
    if not os.path.isfile(Ncontig_return):
        # Without JVM wrapper  
        #pilo = Popen(["java","-Xmx"+args.ram+"G","-jar",Pilon_path+"/pilon-1.23-1.jar","--genome",contig,"--bam",sortedout, "--output",Ncontig,"--outdir",pdir,"--threads",str(args.thread)],stderr=PIPE) 
        # With JVM wrapper  
        pilo = Popen([Pilon_path+"/pilon","-Xmx"+args.ram+"G","--genome",contig,"--bam",sortedout, "--output",Ncontig,"--outdir",pdir,"--threads",str(args.thread)],stderr=PIPE) 
        err=pilo.communicate()
        cmd=" ".join(pilo.args)
        if pilo.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
        logger.debug(cmd)
    Ndid = did + 1
    return (Ncontig_return,Ndid)


def fmlrc(read,correct,did):
    global step
    logger.info(str(step)+". STARTING READS CORRECTION WITH FMLRC")
    time.sleep(2)
    step+=1
    pdir=args.dir+str(did)+"-fmlrc/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    nameR=read.split("/")[-1].split(".")
    nameR=".".join(nameR[:-1])
    read_r=pdir+nameR+"_corrected.fasta"
    nameRea=read.split(".")
    nameRea=".".join(nameRea[:-1])+".fasta"
    if args.retry:
        if not os.path.isfile(pdir+"sr_msbwt.npy") :
            awk=Popen(["awk","NR % 4 == 2",correct],stdout=PIPE)
            sor=Popen(["sort","-T","/datas/Save/Clement/temp"],stdin=awk.stdout,stdout=PIPE)
            tr=Popen(["tr","NT","TN"],stdin=sor.stdout,stdout=PIPE)
            rop=Popen([Ropebwt2_path+"/ropebwt2","-LR"],stdin=tr.stdout,stdout=PIPE)
            tr2=Popen(["tr","NT","TN"],stdin=rop.stdout,stdout=PIPE)
            con=Popen([Fmlrc_path+"/fmlrc-convert","-f",pdir+"sr_msbwt.npy"],stdin=tr2.stdout,stderr=PIPE)
            err=con.communicate()
            cmd=" ".join(awk.args)+" | "+" ".join(sor.args)+" | "+" ".join(tr.args)+" | "+" ".join(rop.args)+" | "+" ".join(tr2.args)+" | "+" ".join(con.args)
            logger.debug(cmd)
            if con.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
        if not os.path.isfile(nameRea) :
            with open(nameRea,"w") as fas:
                sed=Popen(["sed","-n","1~4s/^@/>/p;2~4p",read],stdout=fas)
                sed.wait()
                cmd=" ".join(sed.args)+" "+nameRea
                logger.debug(cmd)
        if not os.path.isfile(read_r) :
            fml=Popen([Fmlrc_path+"/fmlrc","-p",str(args.thread),pdir+"sr_msbwt.npy",nameRea,read_r],stderr=PIPE)
            err=fml.communicate()
            cmd=" ".join(fml.args)
            logger.debug(cmd)
            if fml.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
    else :
        awk=Popen(["awk","NR % 4 == 2",correct],stdout=PIPE)
        sor=Popen(["sort","-T","/datas/Save/Clement/temp"],stdin=awk.stdout,stdout=PIPE)
        tr=Popen(["tr","NT","TN"],stdin=sor.stdout,stdout=PIPE)
        rop=Popen([Ropebwt2_path+"/ropebwt2","-LR"],stdin=tr.stdout,stdout=PIPE)
        tr2=Popen(["tr","NT","TN"],stdin=rop.stdout,stdout=PIPE)
        con=Popen([Fmlrc_path+"/fmlrc-convert","-f",pdir+"sr_msbwt.npy"],stdin=tr2.stdout,stderr=PIPE)
        err=con.communicate()
        cmd=" ".join(awk.args)+" | "+" ".join(sor.args)+" | "+" ".join(tr.args)+" | "+" ".join(rop.args)+" | "+" ".join(tr2.args)+" | "+" ".join(con.args)
        logger.debug(cmd)
        if con.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
        with open(nameRea,"w") as fas:
            sed=Popen(["sed","-n","1~4s/^@/>/p;2~4p",read],stdout=fas)
            sed.wait()
            cmd=" ".join(sed.args)+" "+nameRea
            logger.debug(cmd)
        fml=Popen([Fmlrc_path+"/fmlrc","-p",str(args.thread),pdir+"sr_msbwt.npy",nameRea,read_r],stderr=PIPE)
        err=fml.communicate()
        cmd=" ".join(fml.args)
        logger.debug(cmd)
        if fml.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
    Ndid = did + 1
    return (read_r,Ndid)

def jellyfish(read,did) :
    global step
    pdir=args.dir+str(did)+"-estimate/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    logger.info(str(step)+". STARTING GENOME SIZE ESTIMATION WITH JELLYFISH")
    time.sleep(2)
    step+=1
    if int(args.ram)/args.thread < 1 :
        jelRam = 1
    else :
        jelRam = int(int(args.ram)/args.thread)
    if args.retry :
        if not os.path.isfile(pdir+prefix+"jelly.histo") :
            jel=Popen([Jellyfish_path+"/jellyfish","count","-t",str(args.thread),"-C","-m",str(args.kmer),"-s",str(jelRam)+"G","-o",pdir+prefix+"jelly",read],stderr=PIPE)
            err=jel.communicate()
            cmd=" ".join(jel.args)
            logger.debug(cmd)
            if jel.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
            jelH=Popen([Jellyfish_path+"/jellyfish","histo","-t",str(args.thread),"-o",pdir+prefix+"jelly.histo",pdir+prefix+"jelly"],stderr=PIPE)
            err=jelH.communicate()
            cmd=" ".join(jelH.args)
            logger.debug(cmd)
            if jelH.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
            else :
                try : 
                    os.remove(pdir+prefix+"jelly")
                except Exception as e: 
                    logger.error(e)
                    sys.exit(1)
        if not os.path.isfile(pdir+prefix+"estimate.txt") :
            pass
    else :
        jel=Popen([Jellyfish_path+"/jellyfish","count","-t",str(args.thread),"-C","-m",str(args.kmer),"-s",str(jelRam)+"G","-o",pdir+prefix+"jelly",read],stderr=PIPE)
        err=jel.communicate()
        cmd=" ".join(jel.args)
        logger.debug(cmd)
        if jel.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
        jelH=Popen([Jellyfish_path+"/jellyfish","histo","-t",str(args.thread),"-o",pdir+prefix+"jelly.histo",pdir+prefix+"jelly"],stderr=PIPE)
        err=jelH.communicate()
        cmd=" ".join(jelH.args)
        logger.debug(cmd)
        if jelH.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
        else :
            try : 
                os.remove(pdir+prefix+"jelly")
            except Exception as e: 
                logger.error(e)
                sys.exit(1)
    Ndid = did + 1
    return(Ndid)

def nanoplot(read,did):
    pdir=args.dir+str(did)+"-nanoplot/"
    global step
    logger.info(str(step)+". STARTING READ QUALITY CONTROL WITH NANOPLOT")
    time.sleep(2)
    step+=1
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    if args.retry :
        if not os.path.isfile(pdir+prefix+"-NanoStats.txt") :
            nan=Popen([NanoPlot_path+"/NanoPlot","-t",str(args.thread),"--fastq",read,"-o",pdir,"--N50","--loglength","-p",prefix+"-","--color","cornflowerblue","--store"],stderr=PIPE)
            err=nan.communicate()
            cmd=" ".join(nan.args)
            logger.debug(cmd)
            if nan.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
    else :
        nan=Popen([NanoPlot_path+"/NanoPlot","-t",str(args.thread),"--fastq",read,"-o",pdir,"--N50","--loglength","-p",prefix+"-","--color","cornflowerblue"],stderr=PIPE)
        err=nan.communicate()
        cmd=" ".join(nan.args)
        logger.debug(cmd)
        if nan.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
#        nas=Popen([NanoPlot_path+"/NanoPlot","-t",str(args.thread),"--pickle",pdir+prefix+"-NanoPlot-data.pickle","-o",pdir+"nonlog","-p",prefix+"-","--color","cornflowerblue"],stderr=PIPE)
#        err=nas.communicate()
#        cmd=" ".join(nas.args)
#        logger.debug(cmd)
#        if nas.returncode != 0 :
#            logger.error(err[-1].decode("utf-8"))
#            sys.exit(1)
    Ndid = did + 1
    return (Ndid)


def quast(contig,did):
    pdir=args.dir+str(did)+"-quast/"
    global step
    logger.info(str(step)+". STARTING QUALITY CONTROL WITH QUAST")
    time.sleep(2)
    step+=1
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    if args.retry :
        if not os.path.isfile(pdir+"report.txt") :
            qua=Popen([Quast_path+"/quast",contig,"-o",pdir],stderr=PIPE)
            err=qua.communicate()
            cmd=" ".join(qua.args)
            logger.debug(cmd)
            if qua.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
    else :
        qua=Popen([Quast_path+"/quast",contig,"-o",pdir],stderr=PIPE)
        err=qua.communicate()
        cmd=" ".join(qua.args)
        logger.debug(cmd)
        if qua.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
    Ndid = did + 1
    return (Ndid)


def busco(contig,did):
    pdir=args.dir
    global step
    logger.info(str(step)+". STARTING BUSCO")
    time.sleep(2)
    step+=1
    if not os.path.isdir(pdir):
        os.mkdir(pdir)
    if not os.path.isdir(pdir+str(did)+"-busco/"):
        bus=Popen([Busco_path+"/busco","-l","/home1/datahome/cbellot/augustus_config/eukaryota_odb9/","-i",contig,"--out",timestamp,"--mode","genome","-c",str(args.thread),"-t","./tmp_"+timestamp],stderr=PIPE)
        err=bus.communicate()
        cmd=" ".join(bus.args)
        logger.debug(cmd)
        try :
            shutil.move(os.getcwd()+"/run_"+timestamp+"/short_summary_"+timestamp+".txt",os.getcwd()+"/run_"+timestamp+"/short_summary_busco.txt")
            shutil.move(os.getcwd()+"/run_"+timestamp,pdir+str(did)+"-busco/")
            shutil.move(os.getcwd()+"/tmp_"+timestamp,pdir+str(did)+"-busco/tmp/")
        except : 
            bus=Popen([Busco_path+"/busco","-l","/home1/datahome/cbellot/augustus_config/eukaryota_odb9/","-i",contig,"--out",timestamp,"--mode","genome","-c",str(args.thread)],stdout=PIPE)
            err=bus.communicate()
            print()
            logger.error(err[0].decode("utf-8"))
            sys.exit(1)
    Ndid = did + 1
    return (Ndid) 


def porechop(read,did):
    global step
    logger.info(str(step)+". STARTING TRIMMING WITH PORECHOP")
    time.sleep(2)
    step+=1
    if not args.read :
        logger.error("No read file given (--read).")
        sys.exit(1)
    if "F" not in pat :
        pdir=args.dir+str(did)+"-trim/"
    else :
        pdir=args.dir+str(did)+"-filter/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    nameR=read.split(".")[-2].split("/")[-1]
    read_r=pdir+nameR+"_trimmed.fastq"
    if args.retry :
        if not os.path.isfile(read_r):
            por=Popen([Porechop_path+"/porechop","-t",str(args.thread),"-i",read,"-o",read_r],stderr=PIPE)
            err=por.communicate()
            cmd=" ".join(por.args)
            logger.debug(cmd)
            if por.returncode != 0 :
                logger.error(err[-1].decode("utf-8"))
                sys.exit(1)
    else :
        por=Popen([Porechop_path+"/porechop","-t",str(args.thread),"-i",read,"-o",read_r],stderr=PIPE)
        err=por.communicate()
        cmd=" ".join(por.args)
        logger.debug(cmd)
        if por.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
    if "F" not in pat :
        Ndid = did + 1
    else:
        Ndid = did
    return (read_r,Ndid)


def nanofilt(q,l,read,did):
    global step
    logger.info(str(step)+". STARTING FILTERING WITH NANOFILT")
    time.sleep(2)
    step+=1
    pdir=args.dir+str(did)+"-filter/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    nameR=read.split(".")[-2].split("/")[-1]
    if q is not None and l is not None:
        read_r=pdir+nameR+"_Q"+q+"L"+l+".fastq"
    elif q is not None and l is None:
        read_r=pdir+nameR+"_Q"+q+".fastq"
    elif q is None and l is not None:
        read_r=pdir+nameR+"_L"+l+".fastq"
    else:
        logger.error("ERROR, NanoFilt need at least one threshold (quality or length) to filter the reads.\n")
        sys.exit(1)
    if args.retry :
        if not os.path.isfile(read_r):
            try :
                with open(read,"r") as rea, open(read_r,"w") as reaf :
                    if q is None :
                        nan=Popen([Nanofilt_path+"/NanoFilt","-l",l,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                        err=nan.communicate()
                    elif l is None :
                        nan=Popen([Nanofilt_path+"/NanoFilt","-q",q,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                        err=nan.communicate()
                    else:
                        nan=Popen([NanoFilt_path+"/NanoFilt","-q",q,"-l",l,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                        err=nan.communicate()
                    cmd=" ".join(nan.args)+" "+read+" > "+read_r
                    logger.debug(cmd)
                    if nan.returncode != 0 :
                        logger.error(err[-1].decode("utf-8"))
                        sys.exit(1)
            except Exception as e: 
                logger.error(e)
                sys.exit(1)
    else :
        try :
            with open(read,"r") as rea, open(read_r,"w") as reaf :
                if q is None :
                    nan=Popen([Nanofilt_path+"/NanoFilt","-l",l,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                    err=nan.communicate()
                elif l is None :
                    nan=Popen([Nanofilt_path+"/NanoFilt","-q",q,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                    err=nan.communicate()
                else:
                    nan=Popen([NanoFilt_path+"/NanoFilt","-q",q,"-l",l,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                    err=nan.communicate()
                cmd=" ".join(nan.args)+read+" > "+read_r
                logger.debug(cmd)
                if nan.returncode != 0 :
                    logger.error(err[-1].decode("utf-8"))
                    sys.exit(1)
        except Exception as e: 
            logger.error(e)
            sys.exit(1)
    Ndid = did + 1
    return (read_r,Ndid)

if __name__ == "__main__":

    did = 100
    times=1
    if args.contig :
        contig = args.contig
    else:
        contig = ""
    read = args.read
    correct = args.correct

    if not args.shortread1 :
        sr1 = correct
    else :
        sr1 = args.shortread1

    if not args.shortread2 :
        sr2 = ""
    else :
        if args.shortread1 :
            sr2 = args.shortread2
        else : 
            sr2 = None

    print()
    if not os.path.isdir(args.dir) :
        os.mkdir(args.dir)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    fmt = MyFormatterFile()
    if args.retry :
        logFile = ".".join(logFile.split(".")[:-1])+"_retry.log"
    else :
        logFile = args.dir+"flora.log"
    file_handler = RotatingFileHandler(logFile, "w")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    cfmt = MyFormatterStream()
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(cfmt)
    logger.addHandler(stream_handler)
    if args.retry :
        logger.info("RESTARTING FLORA with "+logFile)
    else : 
        logger.info("STARTING FLORA")
    logger.debug(" ".join(sys.argv))

    for i in pat:
        if i == "E" :
            ret=jellyfish(read,did)
            print()
            did=ret
        elif i == "R" :
            if typ=="fastq" :
                ret=nanoplot(read,did)
                print()
                did=ret
            else :
                logger.info("WARNING, cannot run NanoPlot with file format FASTA. Step ignored.")

        elif i == "F" :
            ret=nanofilt(args.quality,args.length,read,did)
            print()
            read=ret[0]
            did=ret[1]
        elif i == "T":
            ret=porechop(read,did)
            print()
            read=ret[0]
            did=ret[1]
        elif i == "A":
            if not args.assembler :
                sys.exit("ERROR, no assembler given, see --assembler_list for list of available assembler.")
            else:
                if args.assembler=="Flye" or args.assembler=="F":
                    ret=flye(read,did)
                    print()
                    contig=ret[0]
                    did=ret[1]
                elif args.assembler=="Canu" or args.assembler=="C":
                    ret=canu(read,did)
                    print()
                    contig=ret[0]
                    did=ret[1]
                elif args.assembler=="WTDBG2" or args.assembler=="W":
                    ret=wtdbg2(read,did)
                    print()
                    contig=ret[0]
                    did=ret[1]
                elif args.assembler=="Spades" or args.assembler=="S":
                    spades()
                elif args.assembler=="Masurca" or args.assembler=="M":
                    ret=masurca(sr1,sr2,read,did)
                    print()
                    contig=ret[0]
                    did=ret[1]
        elif len(i) == 2 :
            if contig=="" :
                sys.exit("ERROR, no contig file generate or given.")
            if i[0]=="P":
                if args.polisher=="Pilon" or args.polisher=="P" :
                    if i[1]=="s" :
                        ret=pilon("correction","x"+str(times),contig,correct,did)
                        print()
                    else :
                        if read==None:
                            sys.exit("ERROR, no long-read file given.")
                        ret=pilon("polish","x"+str(times),contig,read,did)
                        print()
                elif args.polisher=="Racon" or args.polisher=="R" :
                    if i[1]=="s" :
                        if correct==None:
                            sys.exit("ERROR, no short-read file given.")
                        ret=racon("correction","x"+str(times),contig,correct,did)
                        print()
                    else :
                        if read==None:
                            sys.exit("ERROR, no long-read file given.")
                        ret=racon("polish","x"+str(times),contig,read,did)
                        print()
            times+=1
            contig=ret[0] 
            did=ret[1]
        elif i == "Q":
            ret=quast(contig,did)
            print()
            did=ret
            ret=busco(contig,did)
            print()
            did=ret
        elif i == "C":
            ret=fmlrc(read,correct,did)
            read=ret[0]
            did=ret[1]
            print()

    if "Q" in pat or "R" in pat :
        rmdfile="'"+getReport(args.dir,did,args.fullpath)+"'"
        whi=Popen(["which","Rscript"],stdout=PIPE,stderr=PIPE)
        rsc=whi.communicate()
        if whi.returncode != 0 :
            sys.exit("Rscript not found.")
        else :
            rsc=rsc[0].decode("utf-8") 
        with open(args.dir+"rmark.r","w") as rma :
            rma.write("#! "+rsc)
            rma.write("rmarkdown::render("+rmdfile+")")

        os.chmod(args.dir+"rmark.r", 0o755)
        com=Popen([args.dir+"rmark.r"])
        com.wait()
        os.remove(args.dir+"rmark.r")
    #call('Rscript -e "rmarkdown::render('+report+')"',shell=True)
#    rep=Popen(['Rscript','-e','"rmarkdown::render('+report+')"'])
#    rep.wait()
#    cmd=" ".join(rep.args)
#    print(cmd)
    logger.info("FLORA FINISHED !")

