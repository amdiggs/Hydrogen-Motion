import sys
import os
import re
import linecache as lc
import atom
import numpy as np


#====This is a file that will be good at reading lammps and QE files=============
# ====Steps to making this work well are as follows=========
# 1) get file name and path========
ITEMS=re.compile('ITEM:')
TIMESTEP=re.compile('ITEM: TIMESTEP')
NUM=re.compile('ITEM: NUMBER OF ATOMS')
BB=re.compile('ITEM: BOX BOUNDS')
ATOMS=re.compile('ITEM: ATOMS \w')
Fstart = 0
end = 500
class MD():
    
    def __init__(self):
        self.timestep = None
        self.num = None
        self.box = np.empty([3,2])
        self.atoms = None
        self.done = False



    def set_timestep(self,lines,i):
        self.timestep=int(lines[i + 1])

    def set_num(self,lines,i):
        self.num=int(lines[i + 1])

    def set_bb(self,lines,i):
        for j in range(1,4):
            l=lines[i + j].split()
            self.box[j - 1]=[float(l[0]),float(l[1])]
        print(self.box)   


    def set_atoms(self,lines):
        tmp=np.empty([self.num],dtype=atom.ATOM)
        s=np.diff(self.box,axis=1).reshape(3)
        print(s)
        for line in lines:
           l=line.split()
           x=[int(l[0]),int(l[1]),float(l[2])*s[0],float(l[3])*s[1],float(l[4])*s[2]]
           tmp[x[0]- 1]=atom.ATOM(x)
        self.atoms=tmp
        self.done=True
        return

# I want to read a file
# iterate through each lines
# if ^lines == ITEMS:
# read the next CAPS WORD/S 
# go to switch statment
# if next WORD == TIMESTEP; set timestep
# and so on...
def setup(file,start):
    inst=MD()
    first_line = start
    f=open(file, 'r')
    lines=f.readlines()[start: start + end]
    f.close()
    for i in range(len(lines)):
        if(inst.done):
            if(ITEMS.match(lines[i])):
                break
            else:
                continue
        elif(ITEMS.match(lines[i])):
            set_vals(inst,lines,i)
    atom.ATOM.set_box(inst.box)
    for at in inst.atoms:
        at.get_neighs(inst.atoms)
    global Fstart
    Fstart = first_line + i
    return inst




def set_vals(inst,lines,i):
    if(NUM.match(lines[i])):
        print(lines[i])
        inst.set_num(lines, i)
    elif(BB.match(lines[i])):
        print(lines[i])
        inst.set_bb(lines, i)
    elif(ATOMS.match(lines[i])):
        print(lines[i])
        lst=lines[i+1:i+1+inst.num]
        inst.set_atoms(lst)
    elif(TIMESTEP.match(lines[i])):
        print(lines[i])
        inst.set_timestep(lines, i)
    else:
        print("????")
        print(lines)


### the current issue is reading in multiple time stemps 

def get_neighbors(arr_of_atoms):
#check for tmp_neigh_list.txt
    files=np.asarray(os.listdir("./"))
    if(np.any(np.isin(files,nlist))):
        with open(nlist) as f:
            lines = f.readlines()
            for i in range(len(lines)):
                x=np.array(lines[i].split(),dtype=int)
                arr_of_atoms[i].neighbors = x

    else:
        for at in arr_of_atoms:
            at.get_neighs(arr_of_atoms)

def write_final(moved,ID,infile):
    s=np.diff(atom.ATOM.box,axis=1).reshape(3)
    ids=[]
    for at in moved:
        ids.append(at.id)
    ids=np.asarray(ids)
    at_line="{0} {1} {2:.6f} {3:.6f} {4:.6f}\n"
    ofile='outfiles/final_H_{0}.dat'.format(ID)
    if(len(moved) ==1):
        check = moved[0].distance_moved()
        if(check <= 0.1):
            return None
    
    off=open(ofile, 'w')
    with open(infile) as f:
        for l in f.readlines():
            x=l.split()
            if(ITEMS.match(l)):
                off.write(l)
                continue
            elif(len(x) != 5):
                off.write(l)
                continue
            else:
                rep=np.isin(ids,int(x[0]))
                if(np.any(rep)):
                    at=moved[rep][0]
                    tmp=at.coords/s
                    fline=at_line.format(at.id, at.type, tmp[0], tmp[1], tmp[2])
                    off.write(fline)
                    at.reset()
                    continue
                else:
                    off.write(l)
                    continue
    off.close()
    return ofile

def comp(initial,final):
    print("#####FILE COMP#######")
    s=open(initial, 'r').readlines()
    f=open(final, 'r').readlines()
    for i in range(len(f)):
        if(s[i] != f[i]):
            print(s[i])
            print(f[i])
    print("")
    return 2



if __name__ == "__main__":
    print("???????")



