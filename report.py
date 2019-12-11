import os
import sys
import re
import argparse
import datetime


def parseBusco(flodir,did,fullpath) :   
    dicofile={}
    nbrLine=0
    if flodir[-1] != "/":
        flodir+="/"
    pdir=flodir+str(did)+"-result/"
    if not os.path.isdir(pdir):
        os.makedirs(pdir, exist_ok=True)
    wr="File name\tComplete BUSCO\tComplete and single-copy BUSCO\tComplete and duplicated BUSCO\tFragmented BUSCO\tMissing BUSCO\tTotal BUSCO groups searched\tID\n"
    nbrLine+=1
    businput=[]
    for i in os.walk(flodir):
        for j in i :
            if "short_summary_busco.txt" in j:
                if i[0][-1]!="/" :
                    businput.append(i[0]+"/short_summary_busco.txt")
                else :
                    businput.append(i[0]+"short_summary_busco.txt")
    businput.sort()
    count=0
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
                    if fullpath :
                        wr+="/".join(filename.split("/")[4:6])
                        nbrLine+=1
                    else :
                        wr+=filedir.split("-")[-1].upper()+"_"+str(dicofile[filelogi])
                        nbrLine+=1
                else:
                    res1=re.search("^\s(\d+)\s+\w+",i) # on récupère les chiffres
                    if res1:
                        wr+="\t"+res1.group(1)
                        nbrLine+=1
        count+=1
        wr+="\t"+str(count)+"\n"
        nbrLine+=1
    recapName = pdir+"buscoRecap.csv"
    if nbrLine > 1 :
        with open(recapName,"w") as out:
            out.write(wr)
        return(recapName)
    else :
        return

def parseQuast(flodir,did,fullpath):
    nbrLine=0
    if flodir[-1] != "/":
        flodir+="/"
    pdir=flodir+str(did)+"-result/"
    if not os.path.isdir(pdir):
        os.makedirs(pdir, exist_ok=True)
    wout="File name\t# Contigs >= 500\tLargest contig\tTotal length\tGC (%)\tN50\tID\n"
    nbrLine+=1
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
                    if int(ele[1])==100 :
                        nDidbef = re.search(str(int(ele[1]))+"-\w+",str(i[1]))
                    else:
                        nDidbef = re.search(str(int(ele[1])-1)+"-\w+",str(i[1]))
                    if nDidbef and not fullpath :
                        listName.append(nDidbef.group(0).split("-")[-1])
                    elif nDidbef and fullpath :
                        listName.append(str(i[0]).split("/")[-1]+"/"+nDidbef.group(0))

    dicoName={}
    zipped=list(zip(quainput, listName))
    zipped.sort()
    count=0
    for i,j in zipped :
        with open(i,"r") as report : 
            on = -1
            if j in dicoName.keys() :
                dicoName[j]+=1
            else :
                dicoName[j]=1
            if fullpath :
                wout += j
                nbrLine+=1
            else :
                wout += j.upper()+"_"+str(dicoName[j])
                nbrLine+=1
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
                        nbrLine+=1
                if on == 0 :
                    count+=1
                    wout+="\t"+str(count)+"\n"
                    nbrLine+=1
                    break
    recapName=pdir+"quastRecap.csv"
    if nbrLine > 1 : 
        with open(recapName,"w") as out:
            out.write(wout)
        return(recapName)
    else :
        return

def nanoplotName(flodir,fullpath):
    orde=[]
    for i in os.walk(flodir) :
        if "nanoplot" in i[0] :
            if fullpath :
                orde.append(i[0])
            else :
                val=i[0].split("/")[-2]
                orde.append(i[0].split("/")[-1][0:3]+val)
    orde.sort()
    count=0
    if not fullpath:
        while count < len(orde):
            orde[count]="Nanoplot_"+str(count+1)
            count+=1
#    if not fullpath :
#        x=0
#        while x < len(orde) :
#           orde[x]+=" "+str(x+1) 
#           x+=1
#        x=0
#        while x < len(orde) :
#            orde[x]=orde[x][3:]
#            x+=1
    if len(orde)==1:
        orde[0]=orde[0][:-2]
    return(orde)

def parseNanoplot(flodir,did,fullpath):
    listNanoName=nanoplotName(flodir,fullpath)
    pdir=flodir+str(did)+"-result/"
    if not os.path.isdir(pdir):
        os.makedirs(pdir, exist_ok=True)
    listN=[0]*102
    listUnorder=[]
    listStat=[]
    count=0
    nbrLine=0
    nb=0
    for i in os.walk(flodir):
        if "nanoplot" in i[0] :
            listUnorder.append([i[0],i[2]])
    listUnorder.sort(key = lambda x: x[0])
    for i in listUnorder:
        for j in i[1]:
            if "loglength_kde" in j :
                listN[count+1]=i[0]+"/"+j
                nb+=1
            elif "kde" in j :
                listN[count]=i[0]+"/"+j
                nb+=1
            elif "Weighted_LogTransformed_HistogramReadlength" in j :
                listN[count+5]=i[0]+"/"+j
                nb+=1
            elif "LogTransformed_HistogramReadlength" in j :
                listN[count+3]=i[0]+"/"+j
                nb+=1
            elif "Weighted_HistogramReadlength" in j :
                listN[count+4]=i[0]+"/"+j
                nb+=1
            elif "HistogramReadlength" in j :
                listN[count+2]=i[0]+"/"+j
                nb+=1
            elif "NanoStats" in j :
                listStat.append(i[0]+"/"+j)
        if nb==6 :
            count+=6
            nb=0
    wout="File name\tMean read length\tMean read quality\t# Reads\tN50\tTotal bases\n"
    nbrLine+=1
    countNa=0
    for i in listStat:
        wout+=listNanoName[countNa]+"\t"
        countNa+=1
        with open(i, "r") as nanoStat:
            for j in nanoStat:
                splitJ=j.split(" ")
                if "Mean read length" in j :
                    wout+="".join(splitJ[-1][:-1].split(","))+"\t"
                    nbrLine+=1
                elif "Mean read quality" in j :
                    wout+="".join(splitJ[-1][:-1].split(","))+"\t"
                    #wout+=splitJ[-1][:-1]+"\t"
                    nbrLine+=1
                elif "Number of reads" in j :
                    wout+="".join(splitJ[-1][:-1].split(","))+"\t"
                    nbrLine+=1
                elif "Read length N50" in j :
                    wout+="".join(splitJ[-1][:-1].split(","))+"\t"
                    nbrLine+=1
                elif "Total bases" in j :
                    wout+="".join(splitJ[-1][:-1].split(","))+"\n"
                    nbrLine+=1
    listNano=listN[0:listN.index(0)]
    recapName=pdir+"nanoRecap.csv"
    if nbrLine > 1 :
        with open(recapName,"w") as out:
            out.write(wout)
        return(listNano,listNanoName,recapName)
    else :
        return(listNano,listNanoName,None)

def getReport(gloDir,did,fullpath):
    now = datetime.datetime.now()
    aujour=now.strftime("%d/%m/%y")
    listRep=[]
    recQua=parseQuast(gloDir,did,fullpath)
    recBus=parseBusco(gloDir,did,fullpath)
    parseNano=parseNanoplot(gloDir,did,fullpath)
    listNanoName=parseNano[1]
    listNano=parseNano[0]
    if recQua :
        listRep.append(recQua)
    if recBus :
        listRep.append(recBus)
    pngName=["Kde","Kde_log","ReadLength","ReadLength_log","Weighted_ReadLength","Weighted_ReadLength_log"]
    out="""
---
title: Flora report
output:

    html_document:
        toc: true
        toc_float: true
        df_print: paged
        theme: flatly
---

```{r include = FALSE}
knitr::opts_chunk$set(echo = FALSE)
```
"""
    if listNano :
        out+="""
# Read quality control {.tabset .tabset-fade .tabset-pills}
"""
        if parseNano[2] :
            numline=0
            rowname=""
            rowList=[]
            with open(parseNano[2],"r") as rep :
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
options(DT.options = list(lengthMenu=c(10,20,50), columnDefs = list(list(className ='dt-right', targets =c(1,2,3)),list(className='dt-left',targets=c(0))), language = list(list(search = 'Filter:'),  list(thousand = ','))))

datatable(df, rownames = c({row}), style='bootstrap', class='table-striped table-hover table-bordered', colnames = c({col}), escape = c(FALSE,TRUE))
```
""".format(col=colname, row=rowname)
                out+="""
```{r include = FALSE}
knitr::opts_chunk$set(echo = FALSE)
```
"""
        countPng=0
        for i in listNanoName:
            count=0  
            out+="## "+i+" {.tabset .tabset-pills}\n"
            while count < len(listNano) :
                out+="""
<div id="{idPng}" class="section level3">
### {png}
""".format(png=pngName[count], idPng=pngName[count]+str(countPng))
                out+="![]("+listNano[countPng]+")\n</div>\n"
                countPng+=1
                count+=1
                if countPng % 6 == 0:
                    break
           
    if listRep :
        out+="""
# Assembly quality control {.tabset .tabset-fade .tabset-pills}
"""
        for i in listRep :
            rowname=""
            rowList=[]
            with open(i,"r") as rep :
                if "quast" in i :
                    out+="""
## Quast
    """
                if "busco" in i :
                    out+="""
## Busco
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
    if sys.argv[3]=="True" :
        fullpath=True
    else :
        fullpath=False 
    rmdfile="'"+getReport(sys.argv[1],sys.argv[2],fullpath)+"'"
    whi=Popen(["which","Rscript"],stdout=PIPE,stderr=PIPE)
    rsc=whi.communicate()
    if whi.returncode != 0 :
        sys.exit("Rscript not found.")
    else :
        rsc=rsc[0].decode("utf-8") 
    with open(sys.argv[1]+"/rmark.r","w") as rma :
        rma.write("#! "+rsc)
        rma.write("rmarkdown::render("+rmdfile+")")

    os.chmod(sys.argv[1]+"/rmark.r", 0o755)
    com=Popen([sys.argv[1]+"/rmark.r"])
    com.wait()
    os.remove(sys.argv[1]+"/rmark.r")
