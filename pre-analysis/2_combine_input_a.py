#!/usr/bin/env pypy

import string
import math
import sys 
import os
import numpy as np

#input argument qbox_output.gro file, unit: nm
dt = 10*0.02418885*1.0E-3 #ps - modified

file_num = 50
for i in range(0,8,1):
#list_C = []
#list_O = []
#list_H = []
    List_line_O = []
    List_line_H = []
    List_line_C = []
    total_C = []
    total_O = []
    total_H = []
    for j in range(file_num):
        list_C = []
        list_O = []
        list_H = []
        frame_num_O = 0
        frame_num_H = 0
        frame_num_C = 0
        print("C:",i,"step:",j)
        for line in open("../%d/output_O-%02d.dat"%(j+1, i), 'r'):
            Lwords = line.split()
    
            if Lwords[0] == "Time:":
                #if frame_num_O > 1 or j>0:
                 #   list_O = np.append(list_O, List_line_O, axis=0)
               # if frame_num_O ==0:
                #    pass
                if frame_num_O ==1:
                    list_O = np.append(list_O, List_line_O)
                    list_O = np.array([list_O])
                if frame_num_O > 1:
                    list_O = np.append(list_O, List_line_O, axis=0)

                List_line_O = []
                frame_num_O +=1
            elif len(Lwords)>3:
                Lwords = np.array([float(i) for i in Lwords])
                List_line_O =np.append(List_line_O,Lwords)
                List_line_O=np.array([List_line_O])

        list_O = np.append(list_O, List_line_O, axis=0)
        if j ==0:
            total_O = np.copy(list_O)
        
        
            #print(np.shape(total_O))
        if j > 0:
            while len(total_O[0]) < len(list_O[0]):
                total_O = np.append(total_O, np.zeros((len(total_O),1)), axis=1)
                #print("O:A",j, frame_num_O)
            while len(total_O[0]) > len(list_O[0]):
                list_O = np.append(list_O, np.zeros((len(list_O),1)), axis=1)
                #print("O:B",j, frame_num_O, len(total_O[0]), len(list_O[0]))

            total_O = np.append(total_O, list_O, axis=0)
            
        #print(list_O)
        for line in open("../%d/output_H-%02d.dat"%(j+1, i), 'r'):
            Lwords = line.split()

            if Lwords[0] == "Time:":
                #if frame_num_H ==0:
                 #   pass
                if frame_num_H ==1:
                    list_H = np.append(list_H, List_line_H)
                    list_H = np.array([list_H])
                if frame_num_H >1:
                    list_H = np.append(list_H, List_line_H, axis=0)

                List_line_H = []
                frame_num_H +=1

            elif len(Lwords)>3:
                Lwords = np.array([float(i) for i in Lwords])
                List_line_H=np.append(List_line_H,Lwords)
                List_line_H=np.array([List_line_H])

        list_H = np.append(list_H, List_line_H, axis=0)
        if j ==0:
            total_H = np.copy(list_H)
        
        if j > 0:
            while len(total_H[0]) < len(list_H[0]):
                total_H = np.append(total_H, np.zeros((len(total_H),1)), axis=1)
                #print("H:A",j, frame_num_H)
            while len(total_H[0]) > len(list_H[0]):
                list_H = np.append(list_H, np.zeros((len(list_H),1)), axis=1)
                #print("H:B",j, frame_num_H)
            total_H = np.append(total_H, list_H, axis=0)

        #print(list_O)
        for line in open("../%d/output_C-%02d.dat"%(j+1, i), 'r'):
            Lwords = line.split()

            if Lwords[0] == "Time:":
                #if frame_num_C ==0:
                 #   pass
                if frame_num_C ==1:
                    list_C = np.append(list_C, List_line_C)
                    list_C = np.array([list_C])
                if frame_num_C > 1:
                    list_C = np.append(list_C, List_line_C, axis=0)

                List_line_C = []
                frame_num_C +=1

            elif len(Lwords)>3:
                Lwords = np.array([float(i) for i in Lwords])
                List_line_C=np.append(List_line_C,Lwords)
                List_line_C=np.array([List_line_C])

        list_C = np.append(list_C, List_line_C, axis=0)
        if j ==0:
            total_C = np.copy(list_C)

        if j > 0:
            while len(total_C[0]) < len(list_C[0]):
                total_C = np.append(total_C, np.zeros((len(total_C),1)), axis=1)
                #print("C:A",j,frame_num_C)
            while len(total_C[0]) > len(list_C[0]):
                list_C = np.append(list_C, np.zeros((len(list_C),1)), axis=1)
                #print("C:B",j,frame_num_C)
            total_C = np.append(total_C, list_C, axis=0)



        #print(list_O)
    total_list = np.append(total_O, total_H, axis=1)
    total_list = np.append(total_list, total_C, axis=1)
    for j in range(len(total_list)):
        if(np.isnan(total_list[j,0])):
            total_list[j]=total_list[j-1]

    np.save("outfile-C%02d"%(i),total_list)
