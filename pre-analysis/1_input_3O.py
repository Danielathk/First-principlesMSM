#!/usr/bin/env pypy
'''
input the qbox code output file, and will output each carbon atom's neiboring information
'''
import string
import math
import sys 
import os
import numpy as np

#input argument qbox_output.gro file, unit: nm
dt = 10*0.02418885*1.0E-3 #ps - modified

class atom_type:
    def __init__(self):
    # O, H, C...
        self.xyz = []           # three floats - coordinates
class sort_O_coordinate:         
    def __init__(self, index, Alist, Blist, Clist, a, b, c, index_j, O3nei_list):  #sort_O_coordinate(i, C_list, O_list, H_list, a, b, c, j, carboni.neighbor_Olist)
        self.central_atom=atom_type()     # central atom
        self.neighborO_Clist=[]            # carbon neighbor list
        self.neighborO_Olist=[]            # oxygen neighbor list
        self.neighborO_Hlist=[]            # hydrogen neighbor list
        self.mainframeC=[]
        self.mainframeO_list=[]
        self.index = index
        self.Clist = Alist
        self.Olist = Blist
        self.Hlist = Clist
        self.dimx = a
        self.dimy = b
        self.dimz = c
        self.O_j = index_j
        self.O3_list = O3nei_list[self.O_j][1].xyz
        self.O3_list2 =O3nei_list
    def neiO_Clist(self):
        dist0=distance(self.Clist[self.index].xyz, self.O3_list, self.dimx, self.dimy, self.dimz)
       # print("dist0", dist0)
       # exit()
        for i in self.Clist:
            #if i != self.Clist[self.index]:
               # dist = distance(i.xyz, self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz) 
               # print (dist)
            if (i.xyz[0]-self.O3_list[0])>(self.dimx/2):
                i.xyz[0]=(i.xyz[0]-self.dimx)
            if (i.xyz[0]-self.O3_list[0])<(-self.dimx/2):
                i.xyz[0]=(i.xyz[0]+self.dimx)

            if (i.xyz[1]-self.O3_list[1])>(self.dimy/2):
                i.xyz[1]=(i.xyz[1]-self.dimy)
            if (i.xyz[1]-self.O3_list[1])<(-self.dimy/2):
                i.xyz[1]=(i.xyz[1]+self.dimy)

            if (i.xyz[2]-self.O3_list[2])>(self.dimz/2):
                i.xyz[2]=(i.xyz[2]-self.dimz)
            if (i.xyz[2]-self.O3_list[2])<(-self.dimz/2):
                i.xyz[2]=(i.xyz[2]+self.dimz)
            dist = distance(i.xyz, self.O3_list, self.dimx, self.dimy, self.dimz)
            if dist0==dist:
                self.mainframeC.append([dist, i])
            else:
                self.neighborO_Clist.append([dist, i])


    def neiO_Olist(self):
       # distO0=distance(self.O3_list2[0][1].xyz, self.O3_list, self.dimx, self.dimy, self.dimz)
       # distO1=distance(self.O3_list2[1][1].xyz, self.O3_list, self.dimx, self.dimy, self.dimz)
       # distO2=distance(self.O3_list2[2][1].xyz, self.O3_list, self.dimx, self.dimy, self.dimz)

        #print("d0",distO0, "d1",distO1,"d2",distO2)
        #exit()
        for i in self.Olist:
           # if i != self.Olist[self.index]:
               # dist = distance(i.xyz, self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz) 
               # print (dist)
            if (i.xyz[0]-self.O3_list[0])>(self.dimx/2):
               i.xyz[0]=(i.xyz[0]-self.dimx)
            if (i.xyz[0]-self.O3_list[0])<(-self.dimx/2):
                i.xyz[0]=(i.xyz[0]+self.dimx)

            if (i.xyz[1]-self.O3_list[1])>(self.dimy/2):
                i.xyz[1]=(i.xyz[1]-self.dimy)
            if (i.xyz[1]-self.O3_list[1])<(-self.dimy/2):
                i.xyz[1]=(i.xyz[1]+self.dimy)

            if (i.xyz[2]-self.O3_list[2])>(self.dimz/2):
                i.xyz[2]=(i.xyz[2]-self.dimz)
            if (i.xyz[2]-self.O3_list[2])<(-self.dimz/2):
                i.xyz[2]=(i.xyz[2]+self.dimz)
            dist = distance(i.xyz, self.O3_list, self.dimx, self.dimy, self.dimz)
        #    if (dist==distO0) or (dist==distO1) or (dist==distO2):
        #        self.mainframeO_list.append([dist, i])
        #    else:
            self.neighborO_Olist.append([dist, i])

            #if dist >0.05:
             #   self.neighborO_Clist.append([dist, i])

    def neiO_Hlist(self):
        for i in self.Hlist:
            #if i != self.Clist[self.index]:
               # dist = distance(i.xyz, self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz) 
               # print (dist)
            if (i.xyz[0]-self.O3_list[0])>(self.dimx/2):
                i.xyz[0]=(i.xyz[0]-self.dimx)
            if (i.xyz[0]-self.O3_list[0])<(-self.dimx/2):
                i.xyz[0]=(i.xyz[0]+self.dimx)

            if (i.xyz[1]-self.O3_list[1])>(self.dimy/2):
                i.xyz[1]=(i.xyz[1]-self.dimy)
            if (i.xyz[1]-self.O3_list[1])<(-self.dimy/2):
                i.xyz[1]=(i.xyz[1]+self.dimy)

            if (i.xyz[2]-self.O3_list[2])>(self.dimz/2):
                i.xyz[2]=(i.xyz[2]-self.dimz)
            if (i.xyz[2]-self.O3_list[2])<(-self.dimz/2):
                i.xyz[2]=(i.xyz[2]+self.dimz)
            dist = distance(i.xyz, self.O3_list, self.dimx, self.dimy, self.dimz)

            self.neighborO_Hlist.append([dist, i])


class sort_coordinate:
    def __init__(self, index, Alist, Blist, Clist, a, b, c):
        self.central_atom=atom_type()     # central atom
        self.neighbor_Clist=[]            # carbon neighbor list
        self.neighbor_Olist=[]            # oxygen neighbor list
        self.neighbor_Hlist=[]            # hydrogen neighbor list
        self.index = index
        self.Clist = Alist
        self.Olist = Blist
        self.Hlist = Clist
        self.dimx = a 
        self.dimy = b
        self.dimz = c
    
    def nei_Clist(self):
        for i in self.Clist:
            if i != self.Clist[self.index]:
               # dist = distance(i.xyz, self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz) 
               # print (dist)
                if (i.xyz[0]-self.Clist[self.index].xyz[0])>(self.dimx/2):
                    i.xyz[0]=(i.xyz[0]-self.dimx)
                if (i.xyz[0]-self.Clist[self.index].xyz[0])<(-self.dimx/2):
                    i.xyz[0]=(i.xyz[0]+self.dimx)
                
                if (i.xyz[1]-self.Clist[self.index].xyz[1])>(self.dimy/2):
                    i.xyz[1]=(i.xyz[1]-self.dimy)
                if (i.xyz[1]-self.Clist[self.index].xyz[1])<(-self.dimy/2):
                    i.xyz[1]=(i.xyz[1]+self.dimy)
                    
                if (i.xyz[2]-self.Clist[self.index].xyz[2])>(self.dimz/2):
                    i.xyz[2]=(i.xyz[2]-self.dimz)
                if (i.xyz[2]-self.Clist[self.index].xyz[2])<(-self.dimz/2):
                    i.xyz[2]=(i.xyz[2]+self.dimz)
                dist = distance(i.xyz, self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz)    


                self.neighbor_Clist.append([dist, i])
    def nei_Olist(self):
        for i in self.Olist:
            #dist = distance(i.xyz,self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz)
            if (i.xyz[0]-self.Clist[self.index].xyz[0])>(self.dimx/2):
                i.xyz[0]=(i.xyz[0]-self.dimx)
            if (i.xyz[0]-self.Clist[self.index].xyz[0])<(-self.dimx/2):
                i.xyz[0]=(i.xyz[0]+self.dimx)
                
            if (i.xyz[1]-self.Clist[self.index].xyz[1])>(self.dimy/2):
                i.xyz[1]=(i.xyz[1]-self.dimy)
            if (i.xyz[1]-self.Clist[self.index].xyz[1])<(-self.dimy/2):
                i.xyz[1]=(i.xyz[1]+self.dimy)

            if (i.xyz[2]-self.Clist[self.index].xyz[2])>(self.dimz/2):
                i.xyz[2]=(i.xyz[2]-self.dimz)
            if (i.xyz[2]-self.Clist[self.index].xyz[2])<(-self.dimz/2):
                i.xyz[2]=(i.xyz[2]+self.dimz)
            dist = distance(i.xyz,self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz)


            self.neighbor_Olist.append([dist,i])
    def nei_Hlist(self):
         for i in self.Hlist:
     #        dist = distance(i.xyz,self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz)
    #         print(dist)
             if (i.xyz[0]-self.Clist[self.index].xyz[0])>(self.dimx/2):
                 i.xyz[0]=(i.xyz[0]-self.dimx)
             if (i.xyz[0]-self.Clist[self.index].xyz[0])<(-self.dimx/2):
                 i.xyz[0]=(i.xyz[0]+self.dimx)
                 
             if (i.xyz[1]-self.Clist[self.index].xyz[1])>(self.dimy/2):
                 i.xyz[1]=(i.xyz[1]-self.dimy)
             if (i.xyz[1]-self.Clist[self.index].xyz[1])<(-self.dimy/2):
                 i.xyz[1]=(i.xyz[1]+self.dimy)

             if (i.xyz[2]-self.Clist[self.index].xyz[2])>(self.dimz/2):
                 i.xyz[2]=(i.xyz[2]-self.dimz)
             if (i.xyz[2]-self.Clist[self.index].xyz[2])<(-self.dimz/2):
                 i.xyz[2]=(i.xyz[2]+self.dimz)

             dist = distance(i.xyz,self.Clist[self.index].xyz, self.dimx, self.dimy, self.dimz)
   #          print(dist)
             #print(dist, self.dimx,self.Clist[self.index].xyz, i.xyz)
             self.neighbor_Hlist.append([dist,i])

        
        #return self.neighbor_Clist
    #nei_Clist(self)

#xxx = sort_coordinate(XXXXXXXX)
#xxx.nei_Clist()
#clist = xxx.neighbor_Clist
def distance(O1, O2, a, b, c): 

    if (O1[0]-O2[0] > a/2):
        O1[0]= O1[0]-a
    if (O1[0]-O2[0] < -a/2):
        O1[0]= O1[0]+a

    if (O1[1]-O2[1] > b/2):
        O1[1]= O1[1]-b
    if (O1[1]-O2[1] < -b/2):
        O1[1]= O1[1]+b

    if (O1[2]-O2[2] > c/2):
        O1[2]= O1[2]-c
    if (O1[2]-O2[2] < -c/2):
        O1[2]= O1[2]+c


    xx = abs(O1[0]-O2[0])
    yy = abs(O1[1]-O2[1])
    zz = abs(O1[2]-O2[2])
    return math.sqrt(xx**2 + yy**2 + zz**2)


def rotationMatrix (i, Clist, NeighborOlist,NeighborHlist,a,b,c):
    centralCarbon=np.array(Clist[i].xyz)
    center = np.array(NeighborOlist[0][1].xyz)
    secondnb = np.array(NeighborHlist[0][1].xyz)
    if (centralCarbon[0]-center[0]> a/2):
        centralCarbon[0]=centralCarbon[0]-a
    if (centralCarbon[0]-center[0]< -a/2):
        centralCarbon[0]=centralCarbon[0]+a

    if (centralCarbon[1]-center[1]> b/2):
        centralCarbon[1]=centralCarbon[1]-b
    if (centralCarbon[1]-center[1]< -b/2):
        centralCarbon[1]=centralCarbon[0]+b


    if (centralCarbon[2]-center[2]> c/2):
        centralCarbon[2]=centralCarbon[2]-c
    if (centralCarbon[2]-center[2]< -c/2):
        centralCarbon[2]=centralCarbon[2]+c

#    print("central Oatom:", center)
#    print("centralCarbon:",centralCarbon)
#    print("2ndcloest neighbor:", secondnb)
    
    ria = center-centralCarbon
    rib = center-secondnb
#    print("ria:",ria)
#    print("rib:",rib)
#    print("a,b,c:",a,b,c)
#    exit()
    e_x = ria / np.linalg.norm(ria)
    m_y = rib-np.dot(e_x,rib)*e_x
    m_z = np.cross(ria, rib)

    #e_x = ria / np.linalg.norm(ria)
    e_y = m_y / np.linalg.norm(m_y)
    e_z = m_z / np.linalg.norm(m_z)
    #print(e_x,e_y,e_z,np.dot(e_x,e_y))
    rot_m = np.vstack((e_x, e_y, e_z))
    #print(rot_m)
    return (rot_m.transpose())
#def findneighbors(index,Alist, Blist, Clist, a, b, c):
#    r_cut = min(a/2, b/2, c/2)
    #print ("called:",index)   
    
def transform( rot_M, nei_Clist, nei_Olist, nei_Hlist, a, b, c):
    rij_C = []
    rij_O = []
    rij_H = []
    #print (np.linalg.norm(rot_M))
    r_cut = min(a/2, b/2, c/2)
    for i in range(len(nei_Clist)):
        if nei_Clist[i][0] < r_cut:
            rij = np.array(nei_Olist[0][1].xyz)-np.array(nei_Clist[i][1].xyz)
            rij_tf = np.dot(rij, rot_M)
            #print (rij_tf)
            rij_tf_n = np.linalg.norm(rij_tf)
            rij_temp = np.hstack((rij_tf_n, rij_tf))/(100*rij_tf_n**3)
            rij_C.append(np.array(rij_temp))
          #  print (nei_Clist[i][0], "transform:", rij_C)

    for i in range(len(nei_Olist)):
        if i !=0:
            if nei_Olist[i][0] < r_cut:
                rij = np.array(nei_Olist[0][1].xyz)-np.array(nei_Olist[i][1].xyz)
                rij_tf = np.dot(rij, rot_M)
                #print (rij_tf)
                rij_tf_n = np.linalg.norm(rij_tf)
                rij_temp = np.hstack((rij_tf_n, rij_tf))/(100*rij_tf_n**3)
                rij_O.append(np.array(rij_temp))
       # rij_O[0][0]=0
       # rij_O[0][1]=0
       # rij_O[0][2]=0
       # rij_O[0][3]=0  
          #  print (nei_Olist[i][0], "transform:", rij_O)
  
    for i in range(len(nei_Hlist)):
        if nei_Hlist[i][0] < r_cut:
            rij = np.array(nei_Olist[0][1].xyz)-np.array(nei_Hlist[i][1].xyz)
            rij_tf = np.dot(rij, rot_M)
            #print (rij_tf)
            rij_tf_n = np.linalg.norm(rij_tf)
            rij_temp = np.hstack((rij_tf_n, rij_tf))/(100*rij_tf_n**3)
            rij_H.append(np.array(rij_temp))
          #  print (nei_Hlist[i][0], "transform:", rij_H)
 #   rij_C = np.array(rij_C)
  #  rij_O = np.array(rij_O)
   # rij_H = np.array(rij_H)

    return rij_C, rij_O, rij_H

def distanceH2CO(neibO_H, i, Clist, neibO_Olist,a,b,c):
#(Oa.neighborO_Hlist, i, C_list, carboni.neighbor_Olist)
    Hxyz = np.array(neibO_H[0][1].xyz)
    O1st = np.array(neibO_Olist[0][1].xyz)
    O2nd = np.array(neibO_Olist[1][1].xyz)
    O3rd = np.array(neibO_Olist[2][1].xyz)
    Ccentral = np.array(Clist[i].xyz)
    d1=distance(Ccentral,Hxyz,a,b,c)
    d2=distance(O1st,Hxyz,a,b,c)
    d3=distance(O2nd,Hxyz,a,b,c)
    d4=distance(O3rd,Hxyz,a,b,c)
    if (d2 >= d3) and (d3>= d4):
        ds = np.hstack((0.01/(d1*d1), 0.01/(d4*d4), 0.01/(d3*d3), 0.01/(d2*d2)))
    if (d2 >= d3) and (d3 < d4) and (d4<d2):
        ds = np.hstack((0.01/(d1*d1), 0.01/(d3*d3), 0.01/(d4*d4), 0.01/(d2*d2)))
    if (d2 >= d3) and (d2 <= d4):
        ds = np.hstack((0.01/(d1*d1), 0.01/(d3*d3), 0.01/(d2*d2), 0.01/(d4*d4)))
 
    if (d3 > d2) and (d2>= d4):
        ds = np.hstack((0.01/(d1*d1), 0.01/(d4*d4), 0.01/(d2*d2), 0.01/(d3*d3)))
    if (d3 > d2) and (d3 <= d4):
        ds = np.hstack((0.01/(d1*d1), 0.01/(d2*d2), 0.01/(d3*d3), 0.01/(d4*d4)))
    if (d3 > d2) and (d3 > d4) and (d2 < d4):
        ds = np.hstack((0.01/(d1*d1), 0.01/(d2*d2), 0.01/(d4*d4), 0.01/(d3*d3)))


   # print(ds)
    #exit()
    return ds

def distanceO2OC(neibO_Olist,neibO_Clist, i, Clist, oi, neib_Olist,a,b,c):
#(Oa.neighborO_Olist,Oa.neighborO_Clist, i, C_list, 1, carboni.neighbor_Olist,a,b,c)
    Oxyz = np.array(neib_Olist[oi][1].xyz)
    O1st = np.array(neib_Olist[0][1].xyz)
    O2nd = np.array(neib_Olist[1][1].xyz)
    O3rd = np.array(neib_Olist[2][1].xyz)
    d1=distance(O1st,Oxyz,a,b,c)
    d2=distance(O2nd,Oxyz,a,b,c)
    d3=distance(O3rd,Oxyz,a,b,c)
    if (d1>=d2) and (d2>=d3):
        dO0=d2
        dO1=d1
    if (d1>=d2) and (d3>=d1):
        dO0=d1
        dO1=d3
    if (d1>=d2) and (d2<d3) and (d1>d3):
        dO0=d3
        dO1=d1

    if (d1<d2) and (d2<=d3):
        dO0=d2
        dO1=d3
    if (d1<d2) and (d3<=d1):
        dO0=d1
        dO1=d2
    if (d1<d2) and (d1<d3) and (d2>d3):
        dO0=d3
        dO1=d2

    Ccentral = np.array(Clist[i].xyz)
    dC0=distance(Ccentral,Oxyz,a,b,c)
    dC1=distance(np.array(neibO_Clist[0][1].xyz),Oxyz,a,b,c)
    if dC1 == dC0:
        dC1=distance(np.array(neibO_Clist[1][1].xyz),Oxyz,a,b,c)
   

    ds = np.hstack((0.01/(dO0*dO0), 0.01/(dO1*dO1), 0.01/(dC0*dC0), 0.01/(dC1*dC1)))


    #print(ds)
    #exit()
    return ds





#main()
to_write_O = []
to_write_C = []
to_write_H = []
frame_num = 0
for line in open(sys.argv[1], 'r'):

    Lwords = line.split()

    if Lwords[0] == 'MD,':
        O_list = []
        H_list = []
        C_list = []
        Num_tot = 0

    elif len(Lwords) == 1:
        Num_tot = int(line)     # Total no. of atoms

    # Finding the xyz coordinates of atoms of type O, H, C
    elif len(Lwords) > 3:
        if Lwords[1] == 'O':
            atom = atom_type()
            atom.xyz = line[20:44].split()
            atom.xyz = np.array([float(i) for i in atom.xyz])
            O_list.append(atom)

        elif Lwords[1] == 'H':
            atom = atom_type()
            atom.xyz = line[20:44].split()
            atom.xyz = np.array([float(i) for i in atom.xyz])
            H_list.append(atom)

        elif Lwords[1] == 'C':
            atom = atom_type()
            atom.xyz = line[20:44].split()
            atom.xyz = np.array([float(i) for i in atom.xyz])
            C_list.append(atom)

        elif Lwords[1] == 'Na':
            pass
        else:
            sys.exit("Error in reading!")

    # Setting the cell dimensions
    elif len(Lwords) == 3:
        a = float(Lwords[0])
        b= float(Lwords[1])
        c = float(Lwords[2])
    # fold it back to the original cell
        for i in O_list:
           i.xyz[0] = i.xyz[0]%a
           i.xyz[1] = i.xyz[1]%b
           i.xyz[2] = i.xyz[2]%c
        for i in H_list:
           i.xyz[0] = i.xyz[0]%a
           i.xyz[1] = i.xyz[1]%b
           i.xyz[2] = i.xyz[2]%c
        for i in C_list:
           i.xyz[0] = i.xyz[0]%a
           i.xyz[1] = i.xyz[1]%b
           i.xyz[2] = i.xyz[2]%c

        ''' if frame_num < 1:
            for i in H_list:
                print ("H:",i.xyz)
        '''
        #print ("O, H, C:", len(O_list), len(H_list), len(C_list),frame_num)
        if frame_num ==0:
            output_O=[]
            output_C=[]
            output_H=[]
            for i in range(len(C_list)):
                output_O.append([])
                output_C.append([])
                output_H.append([]) 
              #  output_O.append(open("output_O-%02d.dat"%(i), 'w'))
              #  output_C.append(open("output_C-%02d.dat"%(i), 'w'))
              #  output_H.append(open("output_H-%02d.dat"%(i), 'w'))

        # find Os close to Cs
        #O1nearC1_list = findbonds(C_list, O_list, a, b, c)
        for i in range(len(C_list)):
            carboni = sort_coordinate(i, C_list, O_list, H_list, a, b, c)
            carboni.nei_Clist()
            carboni.nei_Olist()
            carboni.nei_Hlist()
            
            carboni.neighbor_Clist.sort(key = lambda x: float (x[0]))
            carboni.neighbor_Olist.sort(key = lambda x: float (x[0]))
            carboni.neighbor_Hlist.sort(key = lambda x: float (x[0]))


            Oa = sort_O_coordinate(i, C_list, O_list, H_list, a, b, c, 0, carboni.neighbor_Olist)
            Oa.neiO_Clist()
            Oa.neiO_Olist()
            Oa.neiO_Hlist()
            Oa.neighborO_Clist.sort(key = lambda x: float (x[0]))
            Oa.neighborO_Olist.sort(key = lambda x: float (x[0]))
            Oa.neighborO_Hlist.sort(key = lambda x: float (x[0]))
            #Oa.mainframeO_list.sort(key = lambda x: float (x[0]))
            #print("C_origin",Oa.neighborO_Clist)
            #print("O_origin",Oa.neighborO_Olist)
            #print("C_2insert",Oa.mainframeC)
            #print("O_2insert",Oa.mainframeO_list)
            Oa.neighborO_Clist.insert(0, Oa.mainframeC[0])
           # Oa.neighborO_Olist.insert(0, Oa.mainframeO_list[0])
           # Oa.neighborO_Olist.insert(1, Oa.mainframeO_list[1])
           # Oa.neighborO_Olist.insert(2, Oa.mainframeO_list[2])

            #print("C_after",Oa.neighborO_Clist)
            #print("O_after",Oa.neighborO_Olist)
            #exit()


            Ob = sort_O_coordinate(i, C_list, O_list, H_list, a, b, c, 1, carboni.neighbor_Olist)
            Ob.neiO_Clist()
            Ob.neiO_Olist()
            Ob.neiO_Hlist()
            Ob.neighborO_Clist.sort(key = lambda x: float (x[0]))
            Ob.neighborO_Olist.sort(key = lambda x: float (x[0]))
            Ob.neighborO_Hlist.sort(key = lambda x: float (x[0]))
            #Ob.mainframeO_list.sort(key = lambda x: float (x[0]))

            Ob.neighborO_Clist.insert(0, Ob.mainframeC[0])
            #Ob.neighborO_Olist.insert(0, Ob.mainframeO_list[0])
            #Ob.neighborO_Olist.insert(1, Ob.mainframeO_list[1])
            #Ob.neighborO_Olist.insert(2, Ob.mainframeO_list[2])


            Oc = sort_O_coordinate(i, C_list, O_list, H_list, a, b, c, 2, carboni.neighbor_Olist)
            Oc.neiO_Clist()
            Oc.neiO_Olist()
            Oc.neiO_Hlist()
            Oc.neighborO_Clist.sort(key = lambda x: float (x[0]))
            Oc.neighborO_Olist.sort(key = lambda x: float (x[0]))
            Oc.neighborO_Hlist.sort(key = lambda x: float (x[0]))
            #Oc.mainframeO_list.sort(key = lambda x: float (x[0]))

            Oc.neighborO_Clist.insert(0, Oc.mainframeC[0])
            #Oc.neighborO_Olist.insert(0, Oc.mainframeO_list[0])
            #Oc.neighborO_Olist.insert(1, Oc.mainframeO_list[1])
            #Oc.neighborO_Olist.insert(2, Oc.mainframeO_list[2])

         #self.mainframeC=[]
       # self.mainframeO_list=[]

           # for j in range(len(Oa.neighborO_Clist)):
            #     print(frame_num, "C:", i, "neighbor_list", Oa.neighborO_Clist[j][0], Oa.neighborO_Clist[j][1].xyz)
           # for j in range(len(Oa.neighborO_Olist)):
            #     print(frame_num, "O:", i, "neighbor_list", Oa.neighborO_Olist[j][0], Oa.neighborO_Olist[j][1].xyz)

            #for j in range(len(Oa.neighborO_Hlist)):
             #    print(frame_num, "H:", i, "neighbor_list", Oa.neighborO_Hlist[j][0], Oa.neighborO_Hlist[j][1].xyz)

 
            rotMa = rotationMatrix(i, C_list, Oa.neighborO_Olist,Oa.neighborO_Hlist,a,b,c)
            aopC, aopO, aopH = transform(rotMa, Oa.neighborO_Clist, Oa.neighborO_Olist, Oa.neighborO_Hlist, a, b , c)
                     
            rotMb = rotationMatrix(i, C_list, Ob.neighborO_Olist,Ob.neighborO_Hlist,a,b,c)
            bopC, bopO, bopH = transform(rotMb, Ob.neighborO_Clist, Ob.neighborO_Olist, Ob.neighborO_Hlist, a, b , c)
           
            rotMc = rotationMatrix(i, C_list, Oc.neighborO_Olist,Oc.neighborO_Hlist,a,b,c)
            copC, copO, copH = transform(rotMc, Oc.neighborO_Clist, Oc.neighborO_Olist, Oc.neighborO_Hlist, a, b , c)
          
            disHa=distanceH2CO(Oa.neighborO_Hlist, i, C_list, carboni.neighbor_Olist,a,b,c)
            disHb=distanceH2CO(Ob.neighborO_Hlist, i, C_list, carboni.neighbor_Olist,a,b,c)
            disHc=distanceH2CO(Oc.neighborO_Hlist, i, C_list, carboni.neighbor_Olist,a,b,c)
            #distance to oxygens closest to central carbon 



# (Oa.neighborO_Olist,Oa.neighborO_Clist, i, C_list, 1, carboni.neighbor_Olist,a,b,c)
            disOa=distanceO2OC(Oa.neighborO_Olist,Oa.neighborO_Clist, i, C_list, 0, carboni.neighbor_Olist,a,b,c)
            disOb=distanceO2OC(Ob.neighborO_Olist,Ob.neighborO_Clist, i, C_list, 1, carboni.neighbor_Olist,a,b,c)
            disOc=distanceO2OC(Oc.neighborO_Olist,Oc.neighborO_Clist, i, C_list, 2, carboni.neighbor_Olist,a,b,c)

            # print("rota:",rotMa)
           # print("rotc:",rotMb)
           # print("rotc:",rotMc)
            #w,v = np.linalg.eig(rotM)
            #print (w,v)
           # print(rotM)
            #if frame_num ==0 :
            #    for j in range(len(opO)):
             #       print(frame_num, "O:", i, "neighbor_list", 0.1/carboni.neighbor_Olist[j][0], opO[j], carboni.neighbor_Olist[j][1].xyz)
           # if frame_num ==0 :
           #     for j in range(len(carboni.neighbor_Olist)):
            #         print(frame_num, "O:", i, "neighbor_list", carboni.neighbor_Olist[j][0], carboni.neighbor_Olist[j][1].xyz)
           # if frame_num ==0 :
            #    for j in range(len(carboni.neighbor_Hlist)):
             #        print(frame_num, "H:", i, "neighbor_list", carboni.neighbor_Hlist[j][0], carboni.neighbor_Hlist[j][1].xyz)
            #if (np.isnan(opO[0][0])):
             #   opO = output_O[i-1]
             #   opC = output_C[i-1]
             #   opH = output_H[i-1]
             #5O,8H,1C neighbors
            while len(aopO)<6:
                aopO.append([0,0,0,0]) 
            while len(aopH)<8:
                aopH.append([0,0,0,0])
            while len(aopC)<3:
                aopC.append([0,0,0,0])

            while len(bopO)<6:
                bopO.append([0,0,0,0])
            while len(bopH)<8:
                bopH.append([0,0,0,0])
            while len(bopC)<3:
                bopC.append([0,0,0,0])

            while len(copO)<6:
                copO.append([0,0,0,0])
            while len(copH)<8:
                copH.append([0,0,0,0])
            while len(copC)<3:
                copC.append([0,0,0,0])



            neiO_a=aopO[0:6]
            neiH_a=aopH[0:8]
            neiC_a=aopC[0:3]

            neiO_b=bopO[0:6]
            neiH_b=bopH[0:8]
            neiC_b=bopC[0:3]

            neiO_c=copO[0:6]
            neiH_c=copH[0:8]
            neiC_c=copC[0:3]
            
            neiO_a.append(disOa)
            neiO_b.append(disOb)
            neiO_c.append(disOc)

            neiH_a.append(disHa)
            neiH_b.append(disHb)
            neiH_c.append(disHc)

            listO_O=neiO_a+neiO_b+neiO_c
            listO_C=neiC_a+neiC_b+neiC_c
            listO_H=neiH_a+neiH_b+neiH_c
           # 4*3 element are put in H
           # print("O+:",listO_O)
           # print("H+:",listO_H)
           # exit()
           # print("C+:",listO_C)

            output_O[i].append(listO_O)
            output_C[i].append(listO_C)
            output_H[i].append(listO_H) 
	    

         #   if frame_num < 1:
          #       print("Time:", frame_num*dt, aopH)
        frame_num += 1
    
#length_C = []
#length_O = []
#length_H = []
#for i in range(len(C_list)):
 #   for j in range(frame_num):
  #      length_C.append(len(output_C[i][j])) #print("Time:", j*dt, len(output_O[i][j]), len(output_C[i][j]),len(output_H[i][j]))
   #     length_O.append(len(output_O[i][j]))
    #    length_H.append(len(output_H[i][j]))
#for i in range(len(C_list)):
 #   for j in range(frame_num):
  #      print("C:", i, "Time:", j*dt, len(output_O[i][j]), len(output_C[i][j]),len(output_H[i][j]))


for i in range(len(C_list)):
    outputC = open("output_C-%02d.dat"%(i), 'w')
    for j in range(frame_num):
        outputC.write("Time:  %12.8f\n"%(j*dt))
        for k in range(len(output_C[0][0])):
            outputC.write("%12.8f   %12.8f   %12.8f   %12.8f\n"%(output_C[i][j][k][0], output_C[i][j][k][1],output_C[i][j][k][2],output_C[i][j][k][3]))

    outputO = open("output_O-%02d.dat"%(i), 'w')
    for j in range(frame_num):        
        outputO.write("Time:  %12.8f\n"%(j*dt))
        for k in range(len(output_O[0][0])):
            outputO.write("%12.8f   %12.8f   %12.8f   %12.8f\n"%(output_O[i][j][k][0], output_O[i][j][k][1],output_O[i][j][k][2],output_O[i][j][k][3]))

    outputH = open("output_H-%02d.dat"%(i), 'w')
    for j in range(frame_num):        
        outputH.write("Time:  %12.8f\n"%(j*dt))
        for k in range(len(output_H[0][0])):
            outputH.write("%12.8f   %12.8f   %12.8f   %12.8f\n"%(output_H[i][j][k][0], output_H[i][j][k][1],output_H[i][j][k][2],output_H[i][j][k][3]))


