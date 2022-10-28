# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 13:55:15 2018

@author: Administrator
"""
import re
import os
import sys

def calcu_value(value):
    nzz = len(value) - 1
    number = 0
    for i in value[1:]:
        if i == '0': number = number + 1
    #
    per0 = float(number)/nzz
    return per0    

def delete(inputs,outputs):
    of = open(outputs,'a+')
    with open(inputs) as f:
        header = f.readline()
        of.write(header)
        for i in f:
            line = i.replace('\r','').strip('\n').split('\t')
            per0 = calcu_value(line)
            if per0 >= 0.1:
                print str(per0)+'\tdelete!'
            else:
                print str(per0)+'\t not delete!'
                of.write(i)

if __name__ == '__main__':
    inputs = sys.argv[1]
    outputs = sys.argv[2]
    if os.path.exists(inputs):
        if os.path.exists(outputs):
            print outputs+' : file exists!'
            exit()
        else:
            delete(inputs,outputs)
     
