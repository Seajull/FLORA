import os
import sys
import time
import argparse
import re
import logging
from logging.handlers import RotatingFileHandler
from subprocess import call, Popen, PIPE, DEVNULL
from Bio import SeqIO

'''
FLORA stands for (pipeline) For LOng Read Assembly
'''


PARSER = argparse.ArgumentParser()
PARSER.add_argument("-i", "--contig", dest="contig", default=None, help="Input contig files.")
PARSER.add_argument("-r", "--read", dest="read", default=None, help="Input long-read in fastq format. WARNING, you must have the read in fasta format with the exact same name in the same directory than the fastq file.")
PARSER.add_argument("-c", "--correction", dest="correct", default=None, help="Input read used for correction")
PARSER.add_argument("-a", "--assembler", dest="assembler", default=None, choices=["Flye","F","WTDBG2","W","Spades","S"], help="ID of the assembler you want to use (option --assembler_list or -al for the list)")
PARSER.add_argument("-p", "--pattern", dest="pattern", default=None, help="Launch tools in the order of the pattern (ex: FTAQBRcQBRcQB). See --tuto (-u) for a little help and list of letter available")
PARSER.add_argument("-u", "--tuto", dest="tuto", default=None, action="store_true", help="Display a little help and list of letter available for the pattern option")
PARSER.add_argument("-t", "--thread", dest="thread", default=8, type=int, help="Max number of CPU thread available for the tools (default = 8).")
PARSER.add_argument("-m", "--ram", dest="ram", default="64", help="Max giga of RAM available for the tools (default = 64G).")
PARSER.add_argument("-q", "--quality", dest="quality", default=None, help="Quality threshold for filter by NanoFilt (only usefull if Filter step (F) is in pattern).")
PARSER.add_argument("-l", "--length", dest="length", default=None, help="Length threshold for filter by NanoFilt (only usefull if Filter step (F) is in pattern).")
PARSER.add_argument("-e", "--estimate", dest="estimate", default=None, help="Estimate size of the input genome (ex: 80m).")
PARSER.add_argument("-al", "--assembler_list", dest="assembler_list", default=None, action="store_true", help="Display the list of ID of assembler available")
PARSER.add_argument("-o", "--output", dest="output", default="run", help="ID of the run. It will use this ID for output prefix.")
PARSER.add_argument("-d", "--dir", dest="dir", default="./FLORA_OUT/", help="Directory to stock result (default = ./FLORA_OUT/).")

if __name__ == "__main__":

    if len(sys.argv)==1 :
        print()
        PARSER.print_help(sys.stderr)
        print()
        sys.exit(1)
    args = PARSER.parse_args()

    if args.dir[-1] != "/":
        args.dir+="/"

    if args.dir[-1]!="/":
        args.dir+="/"
    if args.dir[0]!="/":
        if args.dir[0]!=".":
            args.dir="./"+args.dir
        elif args.dir[1]!="/":
            args.dir="./"+args.dir[1:]

# Custom formatter
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



    assemblerList = ["Flye","F","WTDBG2","W","Spades","S"]

    global step 
    step=1
    if args.assembler_list: 
        i = 0
        print("\n\tASSEMBLER LIST\n")
        while i < len(assemblerList)-1:
            print('{:.<30s}{:>1s}'.format(assemblerList[i],assemblerList[i+1]))
            i += 2
        print()
        sys.exit()

    patternList = ["Trim","Filter","Correct","Assemble","Quast","Busco","Racon correct","Racon polish","Nanopolish correct","Nanopolish polish","Pilon correct","Pilon polish"]
    letterList = ["T","F","C","A","Q","B","Rc","Rp","Nc","Np","Pc","Pp"]

    if args.tuto:
        i = 0
        print("\n\tTUTORIAL")
        print("\nFLORA need a pattern to known which tools to launch when.\nFor example the pattern AQB launch a Assembly then Quast and then Busco.\n\nYou can chain multiple times the same tool.\nFor example, the pattern FARpQBRcQB launch these step in this order : \n\n\t- Filter \n\t- Assembly \n\t- Racon (polishing with long-read) \n\t- Quast \n\t- Busco \n\t- Racon (correction with short-read) \n\t- Quast \n\t- Busco \n\nSee list below for available option :")
        print()
        while i < len(patternList):
            print('{:.<30s}{:>1s}'.format(patternList[i],letterList[i]))
            i += 1
        print()
        sys.exit()

    if args.output:
        prefix = args.output+"_FLORA"


    if args.correct is None:
        args.correct = args.read

# je sais pas si c'est utile de détecter le format d'entrée des reads mais au cas où
# c'est stocké dans 'typ'
    if args.read:
        ext = args.read.split(".")[-1]
        if ext in ["fa", "fasta"]:
            typ = "fasta"
        elif ext in ["fq", "fastq"]:
            typ = "fastq"
        else:
            sys.exit("ERROR, wrong format for read input (fasta or fastq only).")
    else:
        pass
        #sys.exit("ERROR, need a read file (fasta or fastq).")

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

    if pat.count("A") > 1:
        sys.exit("ERROR, can't do more than one assembly step (A).")

# Comment the following 4 lines if you want to do multiple filter and/or trim step
    if pat.count("T") > 1:
        sys.exit("ERROR, can't do more than one trim step (T).")
    if pat.count("F") > 1:
        sys.exit("ERROR, can't do more than one filter step (F).")

# Not sure about this one, temporarily fix some issue with launching quast, busco, etc before the assembly steps but don't cover all the case (QBA still crash the pipeline, but at least QBPpA not). 
    if "A" in pat[4:]:
        sys.exit("ERROR, assembly step way too late in pattern.")


    gt=round(args.thread/4*3)
    pt=round(args.thread/4)

def flye(read,did):
    global step
    pdir=args.dir+str(did)+"-flye/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    outP=pdir+prefix
    logger.info(str(step)+". STARTING ASSEMBLY USING FLYE")
    time.sleep(3)
    step+=1
    fly=Popen(["flye","--nano-raw",read,"-t",str(args.thread),"--genome-size",args.estimate,"--out-dir",pdir])
    fly.wait()
    cmd=" ".join(fly.args)
    logger.debug(cmd)
#    if fly.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
    contigR = pdir+"scaffolds.fasta"
    Ndid = did + 1
    return (contigR,Ndid)


def wtdbg2(read,did):
    global step
    pdir=args.dir+str(did)+"-wtdbg2/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    outP=pdir+prefix
    logger.info(str(step)+". STARTING ASSEMBLY USING WTDBG2")
    time.sleep(3)
    step+=1
    wtd=Popen(["wtdbg2","-t"+str(args.thread),"-x","ont","-g",args.estimate,"-i",read,"-fo",outP])
    wtd.wait()
    cmd=" ".join(wtd.args)
    logger.debug(cmd) 
#    if wtd.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
    cns=Popen(["wtpoa-cns","-t",str(args.thread),"-i",outP+".ctg.lay.gz","-fo",outP+".ctg.fa"])
    cns.wait()
    cmd=" ".join(cns.args)
    logger.debug(cmd)
#    if cns.returncode != 0 :
#        logger.error(err[-1].decode("utf-8"))
#        sys.exit(1)
    Ndid = did + 1 
    return(outP+".ctg.fa",Ndid)


def spades(sRead1,lRead,sRead2=""):
    print(sRead1)
    print(sRead2)
    print(lRead)
    return


def racon(mode,times,contig,read,did):
    global step
    pdir=args.dir+str(did)+"-racon/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    logger.info(str(step)+". STARTING ALIGNMENT OF CONTIG FILE AND READS")
    time.sleep(3)
    step+=1
    samout = pdir+prefix+"_minimap2_"+times+".sam"
    with open(samout,"w") as alout:
        if mode=="correction":
            mini=Popen(["minimap2","-x","sr","-t",str(gt),"-L","-a",contig, read], stdout=alout,stderr=PIPE)
        else:
            mini=Popen(["minimap2","-t",str(gt),"-L","-a",contig, read], stdout=alout,stderr=PIPE)
        err=mini.communicate()
        cmd=" ".join(mini.args)+" > "+samout
        logger.debug(cmd)
        if mini.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
    logger.info(str(step)+". STARTING RACON "+mode.upper())
    step+=1
    time.sleep(3)
    nameC=contig.split("/")[-1].split(".")
    nameC=".".join(nameC[:-1])
    if mode == "polish" :
        patte=re.compile(".*(_polished_x(\d+)).*")
        mat=patte.match(nameC)
        if mat :
            nameC=nameC.replace(mat.group(1),"_polished_x"+str(int(mat.group(2))+1))
    elif mode == "correction":
        patte=re.compile(".*(_corrected_x(\d+)).*")
        mat=patte.match(nameC)
        if mat :
            nameC=nameC.replace(mat.group(1),"_corrected_x"+str(int(mat.group(2))+1))
    con_r=pdir+nameC+".fasta"
    nameCon=contig.split(".")
    nameCon=".".join(nameCon[:-1])+".fasta"
    with open(con_r,"w") as conout :
        rac = Popen(["racon", "-t",str(args.thread),read,samout,contig], stdout=conout,stderr=PIPE)
    err=rac.communicate()
    print()
    cmd=" ".join(rac.args)+" > "+nameCon
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
    time.sleep(3)
    step+=1
    bamout = pdir+prefix+"_minimap2_"+times+".bam"
    sortedout = pdir+prefix+"_minimap2_sorted_"+times+".bam"
    with open(bamout,"w") as alout:
        if mode=="correction":
            mini=Popen(["minimap2","-x","sr","-t",str(gt),"-L","-a",contig, read], stdout=PIPE)
        else:
            mini=Popen(["minimap2","-t",str(gt),"-L","-a",contig, read], stdout=PIPE)
        samt=Popen(["/home/mbadawi/anaconda3/bin/samtools","view","-@"+str(pt),"-Sb","-"],stdin=mini.stdout,stdout=alout,stderr=PIPE)
        err=samt.communicate()
        cmd=" ".join(mini.args)+" | "+" ".join(samt.args)
        logger.debug(cmd)
        if sams.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
        sams=Popen(["/home/mbadawi/anaconda3/bin/samtools","sort","-@"+str(args.thread),"-o",sortedout,bamout],stderr=PIPE)
        err=sams.communicate()
        cmd=" ".join(sams.args)
        logger.debug(cmd)
        if sams.returncode != 0 :
            logger.error(err[-1].decode("utf-8"))
            sys.exit(1)
        ind=Popen(["/home/mbadawi/anaconda3/bin/samtools", "index","-@"+str(args.thread), sortedout],stderr=PIPE)
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
    time.sleep(3)
    pilo = Popen(["java","-Xmx"+args.ram+"G","-jar","/datas/Save/Clement/soft/pilon/pilon-1.23.jar","--genome",contig,"--bam",sortedout, "--output",Ncontig,"--outdir",pdir,"--threads",str(args.thread)],stderr=PIPE) 
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
    time.sleep(3)
    step+=1
    pdir=args.dir+str(did)+"-fmlrc/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    nameR=read.split("/")[-1].split(".")
    nameR=".".join(nameR[:-1])
    read_r=pdir+nameR+"_corrected.fasta"
    nameRea=read.split(".")
    nameRea=".".join(nameRea[:-1])+".fasta"
    awk=Popen(["awk","NR % 4 == 2",correct],stdout=PIPE)
    sor=Popen(["sort"],stdin=awk.stdout,stdout=PIPE)
    tr=Popen(["tr","NT","TN"],stdin=sor.stdout,stdout=PIPE)
    rop=Popen(["ropebwt2","-LR"],stdin=tr.stdout,stdout=PIPE)
    tr2=Popen(["tr","NT","TN"],stdin=rop.stdout,stdout=PIPE)
    con=Popen(["fmlrc-convert","-f",pdir+"sr_msbwt.npy"],stdin=tr2.stdout,stderr=PIPE)
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
    fml=Popen(["fmlrc","-p",str(args.thread),pdir+"sr_msbwt.npy",nameRea,read_r],stderr=PIPE)
    err=fml.communicate()
    cmd=" ".join(fml.args)
    logger.debug(cmd)
    if fml.returncode != 0 :
        logger.error(err[-1].decode("utf-8"))
        sys.exit(1)
    Ndid = did + 1
    return (read_r,Ndid)


def nanopolish(mode):
    return


def quast(contig,did):
    pdir=args.dir+str(did)+"-quality_control/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    global step
    logger.info(str(step)+". STARTING QUALITY CONTROL WITH QUAST")
    time.sleep(3)
    step+=1
    qua=Popen(["quast.py",contig,"-o",pdir],stderr=PIPE)
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
    if not os.path.isdir(pdir):
        os.mkdir(pdir)
    global step
    logger.info(str(step)+". STARTING BUSCO")
    time.sleep(3)
    step+=1
    bus=Popen(["/datas/Save/Clement/soft/busco/scripts/run_BUSCO.py","-l","/datas/Save/Clement/soft/busco/dataset/eukaryota_odb9/","-i",contig,"--out","busco","--mode","genome","-c",str(args.thread)],stderr=PIPE)
    err=bus.communicate()
    cmd=" ".join(bus.args)
    logger.debug(cmd)
    try :
        os.rename(os.getcwd()+"/run_busco",pdir+str(did)+"-busco/")
        os.rename(os.getcwd()+"/tmp",pdir+str(did)+"-busco/tmp/")
    except : 
        bus=Popen(["/datas/Save/Clement/soft/busco/scripts/run_BUSCO.py","-l","/datas/Save/Clement/soft/busco/dataset/eukaryota_odb9/","-i",contig,"--out","busco","--mode","genome","-c",str(args.thread)],stdout=PIPE)
        err=bus.communicate()
        print()
        logger.error(err[0].decode("utf-8"))
        sys.exit(1)

    Ndid = did + 1
    return (Ndid) 


def porechop(read,did):
    global step
    logger.info(str(step)+". STARTING TRIMMING WITH PORECHOP")
    time.sleep(3)
    step+=1
    if args.read is None :
        logger.error("No read file given (--read).")
        sys.exit(1)
    pdir=args.dir+str(did)+"-filter/"
    if not os.path.isdir(pdir) :
        os.mkdir(pdir)
    nameR=read.split(".")[-2].split("/")[-1]
    read_r=pdir+nameR+"_trimmed.fastq"
    por=Popen(["porechop-runner.py","-t",str(args.thread),"-i",read,"-o",read_r],stderr=PIPE)
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
    time.sleep(3)
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
    try :
        with open(read,"r") as rea, open(read_r,"w") as reaf :
            if q is None :
                nan=Popen(["python3","/datas/Save/Clement/.local/lib/python3.6/site-packages/nanofilt/NanoFilt.py","-l",l,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                err=nan.communicate()
            elif l is None :
                nan=Popen(["python3","/datas/Save/Clement/.local/lib/python3.6/site-packages/nanofilt/NanoFilt.py","-q",q,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
                err=nan.communicate()
            else:
                nan=Popen(["python3","/datas/Save/Clement/.local/lib/python3.6/site-packages/nanofilt/NanoFilt.py","-q",q,"-l",l,"--logfile",pdir+"nanofilt.log"],stdin=rea,stdout=reaf,stderr=PIPE)
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

    #spades("lol","mmh")
    #print("\n")
    #spades("lol1","mmh","lol2")
    did = 100
    times=1
    if args.contig is not None :
        contig = args.contig
    else:
        contig = ""
    read = args.read
    correct = args.correct

    print()
    if not os.path.isdir(args.dir) :
        os.mkdir(args.dir)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    fmt = MyFormatterFile()
    file_handler = RotatingFileHandler(args.dir+"flora.log", "w")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    cfmt = MyFormatterStream()
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(cfmt)
    logger.addHandler(stream_handler)

    logger.info("STARTING FLORA")
    logger.debug(" ".join(sys.argv))

    for i in pat:
        if i == "F" :
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
            if args.assembler is None:
                sys.exit("ERROR, no assembler given, see --assembler_list for list of available assembler.")
            else:
                if args.assembler=="Flye" or args.assembler=="F":
                    ret=flye(read,did)
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
        elif len(i) == 2 :
            if i[0]=="P":
                if i[1]=="c" :
                    ret=pilon("correction","x"+str(times),contig,correct,did)
                    print()
                else :
                    ret=pilon("polish","x"+str(times),contig,read,did)
                    print()
            elif i[0]=="R":
                if i[1]=="c" :
                    ret=racon("correction","x"+str(times),contig,correct,did)
                    print()
                else :
                    ret=racon("polish","x"+str(times),contig,read,did)
                    print()
            else: 
                if i[1]=="c" :
                    ret=nanopolish("correction","x"+str(times),contig,correct,did)
                    print()
                else :
                    ret=nanopolish("polish","x"+str(times),contig,read,did)
                    print()
            times+=1
            contig=ret[0] 
            did=ret[1]
        elif i == "Q":
            ret=quast(contig,did)
            print()
            did=ret
        elif i == "B":
            ret=busco(contig,did)
            print()
            did=ret
        elif i == "C":
            ret=fmlrc(read,correct,did)
            read=ret[0]
            did=ret[1]
            print()
