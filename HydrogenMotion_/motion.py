import readfiles as rf
import atom
import numpy as np
import atom_plot as ap
import matplotlib.pyplot as plt


def region(atom):
    d=5.43
    reg=np.empty([3,2])
    for i in [0,1,2]:
        low, high = atom.coords[i] - d, atom.coords[i] + d
        if(low < atom.box[i][0]):
            low=atom.box[i][0]
        elif(high > atom.box[i][1]):
            high=atom.box[i][1]
        reg[i] = [low,high]
    return reg


def get_env(at):
    i=0
    nlist=[]
    simask=[]
    for n in at.neighbors:
        nlist.append(atom.distance(at,n))
        simask.append(n.type==2)
    nlist=np.asarray(nlist) #array of at neighbor distances
    simask=np.asarray(simask) #Bool array of at silicon neighbors
    bmask=np.logical_and(nlist>1.3,nlist<1.6) #bool array for is at bonded to a Si
    nnmask=np.logical_and(nlist>2.0, nlist<2.8) #bool array for how many neighbors are 2.4->2.8 A away
    sinn=np.logical_and(simask,nnmask) #bool array for are Si neighbors between 2.4->2.8 
    bsi=np.logical_and(simask,bmask) #bool array for is the neighbor 1.4->1.6 a Si atom
    SiBonded=(np.count_nonzero(bsi)==1) #count to make sure only one Si neighbor is 1.4->1.6, if more than 1
    # then H atom is at a bond center.
    count_nn = np.count_nonzero(sinn)
    if(SiBonded and count_nn >= 1):
        nn = at.neighbors[sinn]
        return move_H(at, nn, count_nn)

    else:
        return np.asarray([at])

def move_H(at,nn, cnt):
    if(cnt >= 2):
        for n in nn:
            for oth in nn:
                if(n.id != oth.id):
                    dist=atom.distance(n,oth)
                    if(dist > 2.15 and dist < 2.8):
                        return move_bc(at,n,oth)
    #if no nn are also neighbors then move H to an interstitial site away from the two si atoms.
    return move_inter(at,nn)

def move_bc(at, n1, n2):
    print("BC")
    r=n1.coords - n2.coords
    r_hat=r/np.sqrt(np.sum(r*r))
    mid=(n1.coords + n2.coords)/2.0
    at.coords=mid
    n1.coords= n1.coords + 0.45*r_hat
    n2.coords= n2.coords - 0.45*r_hat
    at.distance_moved()
    n1.distance_moved()
    n2.distance_moved()
    return np.array([at,n1,n2])

def move_inter(at, nn):
    s=np.size(nn)
    if(s >= 2):
        r1= nn[0].coords - at.coords
        r2= nn[1].coords - at.coords
        r_vec= r1 + r2
        r_hat=r_vec/np.sqrt(np.sum(r_vec*r_vec))
        at.coords = at.coords + 3.35*r_hat
        print("2 Neighbors")
        at.distance_moved()
        return np.asarray([at])
    elif(s==1):
        r_vec= nn[0].coords - at.coords
        r_hat=r_vec/np.sqrt(np.sum(r_vec*r_vec))
        at.coords = nn[0].coords + 1.5*r_hat
        print("1 Neighbor")
        at.distance_moved()
        return np.asarray([at])
    else:
        return np.asarray([at])

def get_hydrogen(array_of_atoms):
    mask=np.empty([np.size(ts1.atoms)],dtype=bool)
    i=0
    cont=0
    for at in array_of_atoms:
        if(at.type == 1):
            mask[i]=True
            i+=1
            cont+=1
        elif(at.type == 2):
            mask[i]=False
            i+=1
    return array_of_atoms[mask]

def Migration(at, file):
    print_env(at)
    moved=get_env(at)
    final = rf.write_final(moved, at.id, file)
    if not(final):
        print("Hydrogen atom {0} was not moved")
        return 1
    else:
        rf.comp(file,final)
        return 0

def print_env(at):
    atype=["H","Si"]
    print("")
    txt1="H atom {0} neighbor list".format(at.id)
    txt2="{0} {1} {2:.3f}"
    print(txt1)
    print("Type: ID: Distance:")
    for n in at.neighbors:
        ftext=txt2.format(atype[n.type-1], n.id, atom.distance(at,n))
        print(ftext)
    return 2

def find_inter(at):
    
    return

if __name__ =='__main__':
    FILE="/Users/diggs/Desktop/LAMMPS/AMD_merged/final/final_5.atom"
    ts1= rf.setup(FILE,rf.Fstart)
    hy=get_hydrogen(ts1.atoms)
    i=0
    for h in hy:
        print(i)
        Migration(h,FILE)
        i+=1
















