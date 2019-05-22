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
    with open(cwd+"/config","w") as conf :
        whi=Popen(["which","porechop-runner.py"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Porechop = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","NanoFilt.py"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("NanoFilt = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","ropebwt2"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Ropebwt2 = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","fmlrc"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Fmlrc = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","flye"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Flye = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","wtdbg2"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Wtdbg2 = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","quast.py"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Quast = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","run_BUSCO.py"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Busco = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","racon"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Racon = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","pilon-1.23.jar"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Pilon = "+"/".join(pat.split("/")[:-1])+"\n\n")

        whi=Popen(["which","nanopolish"],stdout=PIPE)
        pat=whi.communicate()
        pat=pat[0].decode("utf-8") 
        conf.write("Nanopolish = "+"/".join(pat.split("/")[:-1])+"\n\n")
        
        print("\n\t"+cwd+"/config file created.\n\n\tPlease check if the path contain in this file are correct.\n")

if __name__ == '__main__' :
    detect()

