# -*- coding: utf-8 -*-
import sys

def setheaderindex(file):
    
    with open(file) as f:
        result = {}
        number = 0
        for i in  f.readline().strip('\n').split('\t'):
            result[i.strip('\r').replace('-','.')] = number
            number = number + 1
        return result

def getcommcol(col1,col2):
    colindex1 = []
    colindex2 = []
    for i in col1.keys():
        if i in col2.keys():
            colindex1.append(col1[i])
            colindex2.append(col2[i])
    return(colindex1,colindex2)
    
def main(m_file,meth_file,out1,out2):
    header1 = setheaderindex(m_file)
    header2 = setheaderindex(meth_file)   
    comm_col = getcommcol(col1=header1,col2=header2)   
    outf1 = open(out1,'a+')
    outf2 = open(out2,'a+')
    with open(m_file) as mf:
        for i in mf:
            ni = i.strip('\r').strip('\n')
            nii = ni.strip('\n').strip('\r')
            context1 = nii.split('\t')
            value = ''        
            for index1 in comm_col[0]:
                value = value +'\t'+ context1[index1]
            line = context1[0]+value+'\n'          
            outf1.write(line)
    print 'over 1'
    kkkkk=1
    with open(meth_file) as mef:
        for j in mef:
            print kkkkk
            kkkkk = kkkkk + 1
            context2 = j.strip('\n').split('\t')
            value = ''
            for index2 in comm_col[1]:
                value = value + '\t' + context2[index2]
            line = context2[0]+'\t'+context2[1]+'\t'+context2[2]+'\t'+context2[3]+'\t'+context2[4]+value+'\n'
            outf2.write(line)
    outf1.close()
    outf2.close()

if __name__ == '__main__':
    if len(sys.argv)<5:
        print "wrong arguments"
        exit(1)
    input1 = sys.argv[1]
    input2 = sys.argv[2]
    output1 = sys.argv[3]
    output2 = sys.argv[4]
    print "input1:"+input1 
    print "input2:"+input2 
    print "output1:"+output1 
    print "output2:"+output2
    main(input1,input2,output1,output2)