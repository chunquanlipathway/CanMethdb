# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 13:55:15 2018

@author: Administrator
"""
import re
import os
import sys
import datetime

def getheader(methfile):
    with open(methfile) as f:
        header = f.readline().replace('\r','').strip('\n').split('\t')
    return header

def getclassindex(header):
    c1index = []
    c2index = []
    number = 0
    for i in header:
        tag = int(re.split('[._]',i.strip('"'))[3])
        if tag > 9:
            c1index.append(number)
        else:
            c2index.append(number)
        number = number + 1
    return (c1index,c2index)

def getclasslocation(sample_class,location):
    if location in sample_class[0]:
        return sample_class[0]
    else:
        return sample_class[1]
def addlist(mlist):
    number = 0.0
    for i in mlist:
        number = number + float(i)
    return number

def deleteNA(methfile,sample_class,out):
    wf = open(out,'a+')
    with open(methfile) as f:
        headerline = f.readline().replace('\r','').strip('\n').split('\t')
        headerkkkk = headerline[0]
        for i in headerline[4:]:
            headerkkkk = headerkkkk +'\t'+i
        wf.write(headerkkkk+'\n')
        linesss_number = 2
        for i in f:
            print str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+' : '+str(linesss_number)
            linesss_number = linesss_number + 1
            line = i.replace('\r','').strip('\n').split('\t')
            value = line[4:]
            line_number = len(value)
            na_number = 0
            for v in  value:
                if v == 'NA':
                    na_number = na_number + 1
                else:
                    pass
            na_v = float(na_number)/line_number
           if na_v == 0:
                new_line = line[0]+'\t'+'\t'.join(value)+'\n'              
                wf.write(new_line)
            elif na_v >=0.1:
                continue
            else:
                na_location = 0
                new_line_value = line[0].strip('"')
                colsss_number = 0
                for na_value in value:
                    colsss_number = colsss_number + 1
                    if na_value == 'NA':
                        m_value = []
                        for cs_vl in getclasslocation(sample_class,na_location):
                             m_value.append(value[cs_vl])  
                        deletena_list = [ m for m in m_value if m!='NA']
                        if len(deletena_list) == 0:
                            break
                        denom = len(deletena_list)
                       m_value = addlist(deletena_list)
                        mean = m_value/denom
                         new_line_value = new_line_value +'\t'+ str(mean)
                    else:
                        new_line_value = new_line_value +'\t'+str(na_value)
                    na_location = na_location + 1
                wf.write(new_line_value+'\n')
    wf.close()

def line100(files):
    wf = open(files+'_100.txt','a+')
    number = 0
    with open(files) as f:
        for i in f:
            if number > 100:
                break
            else:
                wf.write(i)
            number = number + 1
    wf.close()
    
def main(methfile,out):
    header = getheader(methfile)
    sample_class = getclassindex(header[4:])
    deleteNA(methfile,sample_class,out)

if __name__ == '__main__':
    mathfile = sys.argv[1]
    output = sys.argv[2]
    if os.path.exists(output): print 'File! exists!' ;exit()
    main(mathfile,output)
    
else:
    ''''''
    #print '__name__ =>>'+__name__
