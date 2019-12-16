import os
import sys
import argparse
from subprocess import Popen, PIPE

PARSER = argparse.ArgumentParser()
PARSER.add_argument("-f", "--force", dest="force", default=None, action="store_true", help="Force overwrite of config file.")
args = PARSER.parse_args()

def detect():
    cwd = os.path.dirname(os.path.abspath(__file__))
    if not args.force :
        if os.path.isfile(cwd+"/config") :
            sys.exit("Config file already exist, use --force (-f) option to overwrite it.")
    listsoft = ["jellyfish","NanoPlot","porechop","NanoFilt","ropebwt2","fmlrc","flye","canu","wtdbg2","masurca","minimap2","bwa","samtools","pilon","racon","quast","busco"]
    with open(cwd+"/config","w") as conf :
        for i in listsoft :
            whi=Popen(["which",i],stdout=PIPE,stderr=PIPE)
            pat=whi.communicate()
            if whi.returncode == 0 :
                pat=pat[0].decode("utf-8") 
                pat="/".join(pat.split("/")[:-1])
            else :
                if (os.path.isdir("/appli/bioinfo/"+i)):
                    pat="/appli/conda-env/bioinfo/busco-3.0.2/bin"
                else :
                    pat=""
            conf.write(i+" = "+pat+"\n\n")

        print("\n\t"+cwd+"/config file created.\n\n\tPlease check if the path contain in this file are correct.\n")

if __name__ == '__main__' :
    detect()

