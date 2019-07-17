import os
import sys
import re
import datetime


def parseBusco(flodir,did) :   
    dicofile={}
    if flodir[-1] != "/":
        flodir+="/"
    pdir=flodir+str(did)+"-result/"
    if not os.path.isdir(pdir):
        os.makedirs(pdir, exist_ok=True)
    wr="File name\tComplete BUSCO\tComplete and single-copy BUSCO\tComplete and duplicated BUSCO\tFragmented BUSCO\tMissing BUSCO\tTotal BUSCO groups searched\n"
    businput=[]
    for i in os.walk(flodir):
        for j in i :
            if "short_summary_busco.txt" in j:
                if i[0][-1]!="/" :
                    businput.append(i[0]+"/short_summary_busco.txt")
                else :
                    businput.append(i[0]+"short_summary_busco.txt")

    for i in businput :
        with open(i,"r") as summary :
            for i in summary:
                if "Summarized" in i:
                    filename=i.split("file ")[-1][:-1]
                    filedir=filename.split("/")[-2]
                    filelogi=filedir.split("-")[-1]
                    if filelogi not in dicofile.keys() :
                        dicofile[filelogi]=1
                    else:
                        dicofile[filelogi]+=1
                    wr+=filedir.split("-")[-1].upper()+"_"+str(dicofile[filelogi])
                else:
                    res1=re.search("^\s(\d+)\s+\w+",i) # on récupère les chiffres
                    if res1:
                        wr+="\t"+res1.group(1)
        wr+="\n"  
    recapName = pdir+"buscoRecap.csv"
    with open(recapName,"w") as out:
        out.write(wr)
    return(recapName)

def parseQuast(flodir,did):
    if flodir[-1] != "/":
        flodir+="/"
    pdir=flodir+str(did)+"-result/"
    if not os.path.isdir(pdir):
        os.makedirs(pdir, exist_ok=True)
    wout="File name\t# Contigs >= 500\tLargest contig\tTotal length\tGC (%)\tN50\n"
    quainput=[]
    rDid=[]
    listName=[]
    for i in os.walk(flodir):
        if "-quast" in str(i[1]) :
            nDid = re.findall("((\d+)-quast)",str(i[1]))
            if nDid:
                for ele in nDid :
                    for ri in os.walk(i[0]+"/"+ele[0]):
                        if "report.txt" in ri[2] :
                            quainput.append(ri[0]+"/report.txt")
                    nDidbef = re.search(str(int(ele[1])-1)+"-\w+",str(i[1]))
                    if nDidbef :
                        listName.append(nDidbef.group(0).split("-")[-1])
    dicoName={}
    for i,j in zip(quainput, listName) :
        with open(i,"r") as report : 
            on = -1
            if j in dicoName.keys() :
                dicoName[j]+=1
            else :
                dicoName[j]=1
            wout += j.upper()+"_"+str(dicoName[j])
            for i in report :
                if on == -1 :
                    res=re.search("# contigs\s+(\d+)",i)
                    if res: 
                        on = 5
                if on > 0:
                    res=re.search("((\d+\.?\d*)\s*)$",i)
                    on -= 1
                    if res:
                        wout += "\t"+str(res.group(2))
                if on == 0 :
                    wout+="\n"
                    break
    recapName=pdir+"quastRecap.csv"
    with open(recapName,"w") as out:
        out.write(wout)
    return(recapName)

def getReport(gloDir,did):
    now = datetime.datetime.now()
    aujour=now.strftime("%d/%m/%y")
    out="""
---
title: Flora report
output:
    html_document:
        toc: true
        df_print: paged
---

```{r include = FALSE}
knitr::opts_chunk$set(echo = FALSE)
```
## Quality control 

"""
## Quality control {.tabset .tabset-fade}
    listRep=[]
    listRep.append(parseBusco(gloDir,did))
    listRep.append(parseQuast(gloDir,did))
    for i in listRep :
        rowname=""
        rowList=[]
        with open(i,"r") as rep :
            if "quast" in i :
                out+="""
### Quast
"""
            if "busco" in i :
                out+="""
### Busco
"""
            numline=0
            for j in rep :
                if numline < 1 :
                    numline+=1
                    colname=str(j[:-1].split("\t"))[1:-1].replace(", ",",")
                    continue
                rowname+="'<b>"+j.split("\t")[0]+"</b>',"
                rowList.append(str(j[:-1].split("\t")[1:])[1:-1].replace("'","").replace(", ",","))
            rowname=rowname[:-1]

            out+="""
```{r}
library(knitr)
library(DT)
"""
            count=0
            for j in rowList :
                out+="""
r <- list({r})
""".format(r=j)
                if count < 1 :
                    out+="""
df <- data.frame(r,stringsAsFactors = F)
"""
                    count+=1
                else :
                    out+="""
df <- rbind(df,r)
"""
            out+="""
options(DT.options = list(lengthMenu=c(10,20,50), columnDefs = list(list(className ='dt-right', targets =c(1,2,3)),list(className='dt-left',targets=c(0))), language = list(search = 'Filter:')))

datatable(df, rownames = c({row}), style='bootstrap', class='table-striped table-hover table-bordered', colnames = c({col}), escape = c(FALSE,TRUE))
```
""".format(col=colname, row=rowname)
#datatable(df, rownames = c({row}), style='bootstrap', class='table-striped table-hover table-bordered', colnames = c({col}), escape = c(FALSE,TRUE), options = list(order = list(list(0, 'asc'))))
    if gloDir[-1] != "/" :
        gloDir+="/"
    with open (gloDir+"FLORA_report.Rmd","w") as florep :
        florep.write(out)
        print("\n\t---- "+gloDir+"FLORA_report.Rmd created ! ----")
    return(gloDir+"FLORA_report.Rmd")


if __name__=="__main__" :
    from subprocess import call, Popen, PIPE, DEVNULL
    report="'"+getReport(sys.argv[1],sys.argv[2])+"'"
    call('Rscript -e "rmarkdown::render('+report+')"',shell=True)
