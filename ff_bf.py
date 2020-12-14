import pandas as pd 
import random
import numpy

def generate_random_containers(cont,n):
    #ncont no of random containers will be taken from cont for allocation
    ncont=n
    test_containers=[]
    for i in range(ncont):
        rnd=random.randint(0,len(cont)-1)
        while(rnd in test_containers):
            rnd=random.randint(0,len(cont)-1)
        test_containers.append(rnd)
    #print(test_containers)
    return test_containers

def first_fit(cont,vm,pm,test_containers):
    #stats pm=id,spaceavailable,no.of vms
    #pmvm=id,vid,
    #vm=vid,spaceavailable,oid
    #vmcn=vid,cid,
    #cn=cid,oid
    #idcounter=pm,vm,cn

    vmstats=[]
    pmstats=[]
    cnstats=[]
    pmvm=[]
    vmcn=[]

    id_counter=[0,0,0]
    
    for i in test_containers:
        selected_vid=-1
        for vid in range(id_counter[1]):
            ccpu=vmstats[vid][0]-cont[i][0]
            cmem=vmstats[vid][1]-cont[i][1]

            if(ccpu>=0 and cmem>=0):
                selected_vid=vid
                break
        vid=selected_vid
        if(vid!=-1):
            cid=id_counter[2]
            id_counter[2]+=1
            cnstats.append(i)

            vmstats[vid][0]-=cont[i][0]
            vmstats[vid][1]-=cont[i][1]
            vmcn.append([vid,cid])
            vmstats[vid][3]+=1
        else:
                #go over all pm to check availablity of space, create random vm, allocate
                #or create pm, create random vm allocate it here
            selected_vm=random.randint(0,len(vm)-1)
            while(( cont[i][0] > vm[selected_vm][0] ) or ( cont[i][1] > vm[selected_vm][1] ) or ( vm[selected_vm][0] * 1.1 > pm[0] ) or (( vm[selected_vm][1] + 200 ) > pm[1] )):
                selected_vm=random.randint(0,len(vm)-1)

            selected_pid=-1
            for pid in range(id_counter[0]):
            
                if(( pmstats[pid][0] >= vm[selected_vm][0]*1.1 ) and ( pmstats[pid][1] >= vm[selected_vm][1] + 200 )):
                    selected_pid=pid
                    break
                
            pid=selected_pid
            if pid==-1:
                pid=id_counter[0]
                id_counter[0]+=1

                pmstats.append([pm[0],pm[1],0])
                #print(pmstats)
                
            id_counter[1]+=1
            vid=id_counter[1]-1
            vmstats.append([vm[selected_vm][0],vm[selected_vm][1],selected_vm,0])
            
            
            pmstats[pid][0]=pmstats[pid][0]-(vm[selected_vm][0]*1.1)
            pmstats[pid][1]-=vm[selected_vm][1]+200
            pmstats[pid][2]+=1
            pmvm.append([pid,vid])
            
           
            cid=id_counter[2]
            id_counter[2]+=1
           
            cnstats.append(i)
            vmstats[vid][0]-=cont[i][0]
            vmstats[vid][1]-=cont[i][1]
            vmstats[vid][3]+=1
            vmcn.append([vid,cid])
    return pmstats,vmstats,cnstats,pmvm,vmcn

def best_fit(cont,vm,pm,test_containers):
    #stats pm=id,spaceavailable,no.of vms
    #pmvm=id,vid,
    #vm=vid,spaceavailable,oid
    #vmcn=vid,cid,
    #cn=cid,oid
    #idcounter=pm,vm,cn

    vmstats=[]
    pmstats=[]
    cnstats=[]
    pmvm=[]
    vmcn=[]

    id_counter=[0,0,0]
    
    for i in test_containers:
        #searchbf()
        selected_vid=-1
        mincpu=99999
        for vid in range(id_counter[1]):
            ccpu=vmstats[vid][0]-cont[i][0]
            cmem=vmstats[vid][1]-cont[i][1]

            if(ccpu>=0 and cmem>=0 and mincpu>ccpu):
                selected_vid=vid
                mincpu=ccpu
        vid=selected_vid
        if(vid!=-1):
            cid=id_counter[2]
            id_counter[2]+=1
            cnstats.append(i)

            vmstats[vid][0]-=cont[i][0]
            vmstats[vid][1]-=cont[i][1]
            vmcn.append([vid,cid])
            vmstats[vid][3]+=1
        else:
                #go over all pm to check availablity of space, create random vm, allocate
                #or create pm, create random vm allocate it here
            selected_vm=random.randint(0,len(vm)-1)
            while(( cont[i][0] > vm[selected_vm][0] ) or ( cont[i][1] > vm[selected_vm][1] ) or ( vm[selected_vm][0] * 1.1 > pm[0] ) or (( vm[selected_vm][1] + 200 ) > pm[1] )):
                selected_vm=random.randint(0,len(vm)-1)

            selected_pid=-1
            for pid in range(id_counter[0]):
                #if space available
                if(( pmstats[pid][0] >= vm[selected_vm][0]*1.1 ) and ( pmstats[pid][1] >= vm[selected_vm][1] + 200 )):
                    selected_pid=pid
                    break
            pid=selected_pid
            if pid==-1:
                pid=id_counter[0]
                id_counter[0]+=1

                pmstats.append([pm[0],pm[1],0])
                #print(pmstats)
                
            #vmstats,idcounter,vmid=createvm(vmstats,id_counter)
            id_counter[1]+=1
            vid=id_counter[1]-1
            vmstats.append([vm[selected_vm][0],vm[selected_vm][1],selected_vm,0])
            
            
            #allocate_vm(selected_vm,vmid,pid,pmstats,pmvm,vmstats),   also assert if it has the space to do this later
            pmstats[pid][0]=pmstats[pid][0]-(vm[selected_vm][0]*1.1)
            pmstats[pid][1]-=vm[selected_vm][1]+200
            pmstats[pid][2]+=1
            pmvm.append([pid,vid])
            #pmvm[pid][pmstats[pid][2]]=vid

            #allocate_cont(selectedvm,i,vmstats,cnstats),   also assert if it has the space to do this later
            cid=id_counter[2]
            id_counter[2]+=1
            #cnstats[cid]=i
            cnstats.append(i)
            vmstats[vid][0]-=cont[i][0]
            vmstats[vid][1]-=cont[i][1]
            vmstats[vid][3]+=1
            vmcn.append([vid,cid])
    return pmstats,vmstats,cnstats,pmvm,vmcn

#data import 
pm=pd.read_csv('data/PMConfig_small.csv',header=None).to_numpy()
vm=pd.read_csv("data/VMConfig_ten.csv",header=None).to_numpy()
cont=pd.read_csv("data/Container200_ten.csv",header=None).to_numpy()
pmi=[]
for i in range(0,3):
    pmi.append(pm[i][0])
pm=pmi

NO_OF_TEST_CONTAINERS=100

test_containers=generate_random_containers(cont,NO_OF_TEST_CONTAINERS)

pmstats,vmstats,cnstats,pmvm,vmcn=best_fit(cont,vm,pm,test_containers)

print("Best FiT:","\n\n",vmstats)

pmstats,vmstats,cnstats,pmvm,vmcn=first_fit(cont,vm,pm,test_containers)

print("First FiT:","\n\n",vmstats)