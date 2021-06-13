#include<stdio.h>
#include<stdlib.h>
#include<time.h>


#define NO_POP 10                       //controls the no. of solutions to be created in initial population generation
#define NO_CROSSOVER 1                 //controls the no. of times crossover/mutation will be done.
#define CNT_TYPE cont_syn200            //controls the containers array to be used from data1 structure.
#define CNT_SIZE 200                    //same as above, mentions its size.
#define VM_TYPE vm_syn                  //controls the vm array to be used.
#define VM_SIZE 10                      //mentions the vm array size
#define UNPACK_THRSHLD .2               //the threshold below which the pm will be unpacked in crossover/mutation.
#define PROB_PM_TO_ADD 5                //the probability for which best pm is chosen in crossover.
#define PROB_GENE_BIN_TRN 8             //the probability for which best gene is returned in binary tournament
#define PROB_CROSSOVER 10
#define MUTATION_VM_REM 5
#define MUTATION_VM_PROB 8
#define MUTATION_VM_THRSHLD .2

/*
    heuristics:

    1:check if best_pm has space,if yes, use,dbl replacement, complete fit, better fit, best fit in the respective order.
    2: if not do check with the next best_pm.
    3: if none has space, create new.
*/

struct data{
    float cont_syn200[200][2];
    float cont_rea200[200][2];
    float cont_syn500[500][2];
    float cont_rea500[500][2];
    float cont_syn1000[1000][2];
    float cont_rea1000[1000][2];
    float cont_syn1500[1500][2];
    float cont_rea1500[1500][2];
    float vm_syn[10][2];
    float vm_rea[20][2];
    float pm[3];
    int test_containers[100];
    int containers_sample[200][2];
}data1;
struct cn{
    int ogi;//og index
    struct cn *next;
    float cpu;
    float mem;
    float norm_res;

};
struct vm{
    int ogi;
    int ncn;
    float cpu;
    float mem;
    float utlz;
    float mem_utlz;
    struct vm *next;
    struct cn *hcn;
};
struct pm{
    int nvm;
    float cpu;
    float mem;
    float utlz;
    float mem_utlz;
    struct pm *next;
    struct vm *hvm;
};
struct gene{
    float utlz;
    int npm;
    struct pm *hpm;
    int test_containers[CNT_SIZE];
    struct gene *next;
};
struct population{
    struct gene *gene;
    int n;
}*head_population;


typedef struct population* popptr;
typedef struct gene* geneptr;
typedef struct pm* pmptr;
typedef struct vm* vmptr;
typedef struct cn* cnptr;



void readCsv();//reads all the data into data1 structure
void createPop();//creates initial population of solutions
void generate_random_containers(geneptr t1);//creates random iteration of containers to be allocated and stores them in t1->heuristics.
void first_fit(geneptr t1);//allocates containers present in t1->test_containers, acc to first fit heuristics in t1
void printstats(geneptr t1);//prints the stats of t1
vmptr ff_vm(geneptr t1,float cpureq,float memreq);//sub function of first_fit, it finds a vm that fits the cpu and mem req of container.
pmptr ff_pm(geneptr t1,float vmCpu,float vmMem);//sub function of first_fit, it finds a pm that fits cpu and mem req of vm.
pmptr create_pm(geneptr t1);//creates a new pm in t1
vmptr create_vm(int vmi);//creates a new vm
cnptr create_cn(int cni);//creates a new container
void addVmInPm(vmptr vm,pmptr pm);//adds vm into pm
void addCnInVm(cnptr cn,vmptr vm);//adds cn into vm
void fitness(geneptr t1);//calculates the fitness of t1
void search_best_gene(popptr p1);//prints the best gene acc to fitness
geneptr crossover(popptr p1);
geneptr binary_trnment(popptr p1);
pmptr search_best_pm(geneptr t1);
void copy_pm(pmptr pm1,pmptr pm2);
void remove_pm(pmptr pm,geneptr t1);
void removeContainer(int ogi,geneptr t1);
void cnInGene(geneptr t1,int* cn_present);
void count_containers(geneptr t1);
void integrity_all_containers(geneptr t1);
void heuristics(geneptr t1,int cn_absent[CNT_SIZE],int absent_count);
void integrity_duplicate_containers(geneptr t1);
void remove_gene(popptr p1,geneptr t1);
void add_gene(popptr p1,geneptr t1);
geneptr create_gene();
int bf_vm_new(int cn_ogi);
int ff_vm_new(int cn_ogi);
void unpack_pm(geneptr t1);
void merge_vm(geneptr t1);
void remove_vm(pmptr pm,vmptr vm);
void merge2vm(vmptr vm1,vmptr vm2,pmptr pm);
void fitness1(geneptr t1);
void sortPmByUtlz(geneptr t1,int cpumem);
pmptr findprevpm(geneptr t1,pmptr pm);
void switchpm(pmptr pm1,pmptr pm2,geneptr t1);
void crossMut(popptr p1);
int dbl_rep_heuristics(geneptr t1,pmptr pm1,int cn_ogi,int* cn_absent,int* absent_count);
int complete_fill_heuristics(pmptr pm1,int cn_ogi);
int better_fill_heuristics(geneptr t1,pmptr pm1,int cn_ogi,int *cn_absent,int *absent_count);
int best_fill_heuristics(pmptr pm1,int cn_ogi);
int first_fill_heuristics(pmptr pm1, int cn_ogi);
geneptr crossover2(geneptr parent1,geneptr parent2,int cpumem);
void heuristics2(geneptr t1,int cn_absent[CNT_SIZE],int absent_count);
void sortContainersByRes(geneptr t1);
void calculateContRes(geneptr t1);
cnptr findprevcn(cnptr c1,vmptr v1);
int rep_small_two_heuristics(geneptr t1,vmptr v1,int cn_absent[CNT_SIZE],int absent_count);
void ff_rc_heuristics(cnptr cn,geneptr t1);

int main (){

    srand(time(0));
    readCsv();
    printf("\nread done");
    createPop();;
    for(int i=0;i<NO_CROSSOVER;i++){
        crossMut(head_population);
    }
    search_best_gene(head_population);
    printf("\nend\n");
}


/*void unpack_pm(geneptr t1){
    int absentcn[CNT_SIZE];
    int absent_count=0;
    sortPmByUtlz(t1);
    pmptr lowest_pm;
    for(pmptr pm1=t1->hpm;pm1->next!=NULL;pm1=pm1->next){
        lowest_pm=pm1;
    }


    for(vmptr vm1=lowest_pm->hvm;vm1!=NULL;vm1=vm1->next){
        for(cnptr cn1=vm1->hcn;cn1!=NULL;cn1=cn1->next){
            absentcn[absent_count]=cn1->ogi;
            absent_count++;
        }
    }
    remove_pm(lowest_pm,t1);

    heuristics(t1,absentcn,absent_count);
}*/

void unpack_pm(geneptr t1){
    int absentcn[CNT_SIZE];

    int absent_count=0;
    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        if(pm1->utlz>UNPACK_THRSHLD){
            continue;
        }
        for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            for(cnptr cn1=vm1->hcn;cn1!=NULL;cn1=cn1->next){
                absentcn[absent_count]=cn1->ogi;
                absent_count++;
            }
        }
        remove_pm(pm1,t1);
    }
    heuristics(t1,absentcn,absent_count);
}

geneptr create_gene(){
    geneptr t2=(geneptr)malloc(sizeof(struct gene));
    t2->npm=0;
    t2->hpm=NULL;
    t2->utlz=0;
    t2->next=NULL;
    return(t2);

}

void remove_gene(popptr p1,geneptr t1){
    if(p1->n<=0){
        printf("\nERROR:NO GENE PRESENT TO REMOVE");
        exit(0);
    }
    //only 1 present
    if(p1->n==1){
        p1->gene=NULL;
        p1->n=0;
        return;
    }
    //first remove
    geneptr t2=p1->gene;
    if(t2==t1){
        p1->gene=t2->next;
        p1->n--;
        return;
    }

    for(geneptr g1=p1->gene;g1!=NULL;g1=g1->next){
        if(g1->next==t1){
            g1->next=t1->next;
            p1->n--;
            return;
        }
    }
}

void add_gene(popptr p1,geneptr t1){
    if(p1->n==0){
        p1->n++;
        p1->gene=t1;
        return;
    }
    geneptr g1=p1->gene;
    while(g1->next!=NULL){
        g1=g1->next;
    }
    g1->next=t1;
    t1->next=NULL;
    p1->n++;
}

void integrity_duplicate_containers(geneptr t1){
    int cn_present[CNT_SIZE];
    for(int i=0;i<CNT_SIZE;i++){
        cn_present[i]=0;
    }
    pmptr pm1=t1->hpm;
    while(pm1!=NULL){
        vmptr vm1=pm1->hvm;
        while(vm1!=NULL){
            cnptr cn1=vm1->hcn;
            while(cn1!=NULL){
                if(cn_present[cn1->ogi]==1){
                    printf("\n\n\tIntegrity error:DUPLICATE CONTAINERS PRESENT:%d",cn1->ogi);
                    exit(0);
                }
                cn_present[cn1->ogi]=1;
                cn1=cn1->next;
            }
            vm1=vm1->next;
        }
        pm1=pm1->next;
    }
}

void integrity_all_containers(geneptr t1){
    int cn_present[CNT_SIZE];
    int count=0;
    for(int i=0;i<CNT_SIZE;i++){
        cn_present[i]=0;
    }
    pmptr pm1=t1->hpm;
    while(pm1!=NULL){
        vmptr vm1=pm1->hvm;
        while(vm1!=NULL){
            cnptr cn1=vm1->hcn;
            while(cn1!=NULL){
                cn_present[cn1->ogi]=1;

                cn1=cn1->next;
            }
            count+=vm1->ncn;
            vm1=vm1->next;
        }
        pm1=pm1->next;
    }
    if(count!=CNT_SIZE){
        printf("\n\n\tIntegrity error:UNEQUAL NO OF CONTAINERS");
        exit(0);
    }
    for(int i=0;i<CNT_SIZE;i++){
        if(cn_present[i]==0){
            printf("\n\n\tIntegrity error:NOT ALL CONTAINERS PRESENT");
            exit(0);
        }
    }
}

void removeContainer(int ogi,geneptr t1){
    pmptr pm1=t1->hpm;
    //printf("\nhere 1");
    while(pm1!=NULL){
        //printf("\nhere 2");
        vmptr vm1=pm1->hvm;
        while(vm1!=NULL){
            //printf("\nhere 3");
            cnptr cn1=vm1->hcn;
            cnptr cn2;
            //if vm is empty
            if(cn1==NULL){
                vm1=vm1->next;
                continue;
            }
            //if cn is the first container
            if(cn1->ogi==ogi){
                vm1->hcn=cn1->next;
                vm1->cpu+=cn1->cpu;
                vm1->mem+=cn1->mem;
                vm1->ncn--;
                return;
            }
            //else
            while(cn1!=NULL){
                if(cn1->ogi==ogi){
                    vm1->cpu+=cn1->cpu;
                    vm1->mem+=cn1->mem;
                    cn2->next=cn1->next;
                    vm1->ncn--;
                    return;
                }
                cn2=cn1;
                cn1=cn1->next;
            }
            vm1=vm1->next;
        }
        pm1=pm1->next;
    }
}

void remove_pm(pmptr pm,geneptr t1){
    if(t1->hpm==NULL){
        printf("\n PM already empty");
        exit(0);
    }
    pmptr pm1=t1->hpm;
    t1->npm--;
    if(pm1==pm){
        t1->hpm=pm1->next;
        return;
    }
    //printf("\n PM 3");

    while((pm1->next!=pm)){
        pm1=pm1->next;
    }

    pm1->next=pm->next;

    return;

}

void copy_pm(pmptr pm1,pmptr pm2){
    vmptr vm1=pm1->hvm;
    while(vm1!=NULL){
        vmptr vm2=create_vm(vm1->ogi);
        cnptr cn1=vm1->hcn;
        while(cn1!=NULL){
            cnptr cn2=create_cn(cn1->ogi);
            addCnInVm(cn2,vm2);
            cn1=cn1->next;
        }
        addVmInPm(vm2,pm2);
        vm1=vm1->next;
    }
}

void cnInGene(geneptr t1,int* cn_present){
    for(int i=0;i<CNT_SIZE;i++){
        cn_present[i]=0;
    }
    pmptr pm1=t1->hpm;
    while(pm1!=NULL){
        vmptr vm1=pm1->hvm;
        while(vm1!=NULL){
            cnptr cn1=vm1->hcn;
            while(cn1!=NULL){
                if(cn_present[cn1->ogi]==1){
                    printf("\n\n\tIntegrity error:DUPLICATE CONTAINERS PRESENT");
                    exit(0);
                }
                cn_present[cn1->ogi]=1;
                cn1=cn1->next;
            }
            vm1=vm1->next;
        }
        pm1=pm1->next;
    }

}

void count_containers(geneptr t1){
    pmptr pm1=t1->hpm;
    int count=0;
    while(pm1!=NULL){
        vmptr vm1=pm1->hvm;
        while(vm1!=NULL){
            count+=vm1->ncn;
            vm1=vm1->next;
        }
        pm1=pm1->next;
    }
    printf("\nno of containers:%d",count);
}

void merge2vm(vmptr vm1,vmptr vm2,pmptr pm1){
    //merge vm1 into vm2
    //printf("\nmerge2vm 1");
    for(cnptr cn1=vm1->hcn;cn1!=NULL;cn1=cn1->next){
       // printf("\nmerge2vm 2");
        cnptr newCn=create_cn(cn1->ogi);
        addCnInVm(newCn,vm2);
    }
    //printf("\nmerge2vm 3");
    remove_vm(pm1,vm1);
    //printf("\nmerge2vm 4");
}

void remove_vm(pmptr pm,vmptr vm){
    //empty case
    //printf("\nremove_vm 1");
    if((pm->nvm==0)||(pm->hvm==NULL)){
        printf("\nERROR:NO VM FOUND TO REMOVE");
        exit(0);
    }
    //first value case
    //printf("\nremove_vm 2");
    vmptr vm2=pm->hvm;
    float vm_cpu=data1.VM_TYPE[vm->ogi][0];
    float vm_mem=data1.VM_TYPE[vm->ogi][1];
    pm->cpu+=(vm_cpu*1.1);
    pm->mem+=(vm_mem+200);
    //printf("\nremove_vm 3");

    if(vm2==vm){
        //printf("\nremove_vm 4");
        pm->hvm=vm2->next;
        pm->nvm--;
        return;
    }
    //any value
   // printf("\nremove_vm 5");
    while(vm2->next!=vm){
        vm2=vm2->next;
    }
    //printf("\nremove_vm 6");
    if(vm->next==NULL){
        vm2->next=NULL;
    }
    else
        vm2->next=vm->next;
    pm->nvm--;
}

void merge_vm(geneptr t1){
    int count=0;
    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        //int dc=1;
        //vmptr vm1;
        //empty case
        if(pm1->hvm==NULL){
            continue;
        }
        for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            //dc=0;
            //vm1=pm1->hvm;
            float cpu_occ=data1.VM_TYPE[vm1->ogi][0] - vm1->cpu;
            float mem_occ=data1.VM_TYPE[vm1->ogi][1] - vm1->mem;

            for(vmptr vm2=pm1->hvm;vm2!=NULL;vm2=vm2->next){
                if(vm2==vm1){
                    continue;
                }
                if((cpu_occ>vm2->cpu)||(mem_occ>vm2->mem)){
                    continue;
                }
                merge2vm(vm1,vm2,pm1);
                //dc=1;
                count++;
                break;
            }
        }
    }
    printf("\n%d vm merged",count);
    fitness(t1);
}
/*
void merge_vm(geneptr t1){
    int count=0;
    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        //int dc=1;
        //vmptr vm1;
        //empty case
        if(pm1->hvm==NULL){
            continue;
        }
        for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            //dc=0;
            //vm1=pm1->hvm;
            float cpu_occ=data1.VM_TYPE[vm1->ogi][0] - vm1->cpu;
            float mem_occ=data1.VM_TYPE[vm1->ogi][1] - vm1->mem;
            int thisdone=0;
            for(pmptr pm2=t1->hpm;pm2!=NULL;pm2=pm2->next){
                if(thisdone)
                    break;
                for(vmptr vm2=pm2->hvm;vm2!=NULL;vm2=vm2->next){
                   // printf("\nmerge 2");
                    if(vm2==vm1){
                        continue;
                    }
                    if((cpu_occ>vm2->cpu)||(mem_occ>vm2->mem)){
                        continue;
                    }
                    //printstats(t1);
                    //printf("\nmerge 3 %f, %f\t%f, %f\t%d,%d",cpu_occ,vm2->cpu,mem_occ,vm2->mem,pm1->nvm,pm2->nvm);
                    merge2vm(vm1,vm2,pm1);
                    //printf("\nmerge 4");
                    //dc=1;
                    count++;
                    thisdone=1;
                    break;
                }
            }
        }
    }
    printf("\n%d vm merged",count);
    fitness(t1);
}
*/
pmptr findprevpm(geneptr t1,pmptr pm){
    pmptr prev=t1->hpm;
    if(pm==t1->hpm){
        printf("\nERROR:Cant find prev of head pm");
        exit(0);
    }
    while((prev->next!=pm)&&(prev!=NULL)){
        prev=prev->next;
    }
    if(prev==NULL){
        printf("\nERROR:Cant find pm in the given gene");
        exit(0);
    }
    return(prev);

}

void switchpm(pmptr pm1,pmptr pm2,geneptr t1){
    if((pm1==NULL)||(pm2==NULL)){
        printf("\nERROR:Cant switch NULL");
        exit(0);
    }
    if(pm1==pm2){
        return;
    }
    if(pm1->next==pm2){
        if(pm1==t1->hpm){
            //if pm2 is last
            //printf("\nswitch 2");
            t1->hpm=pm2;

            if(pm2->next==NULL){
                  //  printf("\nswitch 3");
                    pm2->next=pm1;
                    pm1->next=NULL;
                    return;
            }
            //printf("\nswitch 4");
            pmptr pm2next=pm2->next;
            pm2->next=pm1;
            pm1->next=pm2next;
            return;

        }

        //pm2 is last but pm1 isnt first
        if(pm2->next==NULL){ //changed pm1 to pm2 here
            //printf("\nswitch 6");
            pmptr prevpm1=findprevpm(t1,pm1);

            prevpm1->next=pm2;
            pm1->next=NULL;
            pm2->next=pm1;
            return;
        }
        //if neither

        pmptr prevpm1=findprevpm(t1,pm1);
        pmptr nextpm2=pm2->next;

        prevpm1->next=pm2;
        pm2->next=pm1;
        pm1->next=nextpm2;
        return;
    }
    if(pm2->next==pm1){
        switchpm(pm2,pm1,t1);
        return;
    }
    //pm1 is first
    //printf("\nswitch 1");
    if(pm1==t1->hpm){
        //if pm2 is last
        //printf("\nswitch 2");
        pmptr pm2prev=findprevpm(t1,pm2);
        t1->hpm=pm2;

        if(pm2->next==NULL){
                //printf("\nswitch 3");
                pm2->next=pm1->next;
                pm2prev->next=pm1;
                pm1->next=NULL;
                return;
        }
        //printf("\nswitch 4");
        pmptr pm2next=pm2->next;
        pm2->next=pm1->next;
        pm2prev->next=pm1;
        pm1->next=pm2next;
        return;

    }
    //pm2 is first
    if(pm2==t1->hpm){
        //printf("\nswitch 5");
        switchpm(pm2,pm1,t1);
        return;
    }

    //pm1 is last
    if(pm1->next==NULL){
       // printf("\nswitch 6");
        pmptr prevpm1=findprevpm(t1,pm1);
        pmptr prevpm2=findprevpm(t1,pm2);
        pmptr nextpm2=pm2->next;

        prevpm1->next=pm2;
        prevpm2->next=pm1;
        pm1->next=nextpm2;
        pm2->next=NULL;
        return;
    }
    //pm2 is last
    if(pm2->next==NULL){
        //printf("\nswitch 7");
        switchpm(pm2,pm1,t1);
        return;
    }
   // printf("\nswitch 8");
    pmptr prevpm1=findprevpm(t1,pm1);
    pmptr prevpm2=findprevpm(t1,pm2);
    pmptr nextpm1=pm1->next;
    pmptr nextpm2=pm2->next;

    prevpm1->next=pm2;
    prevpm2->next=pm1;
    pm1->next=pm2;
    pm2->next=pm1;
    return;
}

void sortPmByUtlz(geneptr t1, int cpumem){

    //can break at the final case when next=NULL
    //printf("\ninRep 1");
    fitness1(t1);
    //printf("\ninRep 2");
    pmptr pm1=t1->hpm;
    if(cpumem==0){
        while(pm1->next!=NULL){
            pmptr temp=pm1;
            //printf("\ninRep 2.1");
            int i=0;
            for(pmptr pm2=pm1->next;pm2!=NULL;pm2=pm2->next){
                i++;
                //printf("\ninRep 2.11 %d",i);
                if(temp->mem_utlz<pm2->mem_utlz){
                    temp=pm2;
                }
            }
            if(temp==pm1){
                pm1=pm1->next;
                continue;
            }
            pmptr pmnext=pm1->next;

            switchpm(pm1,temp,t1);

            pm1=pmnext;
        }
    }
    else{
        while(pm1->next!=NULL){
            pmptr temp=pm1;
            //printf("\ninRep 2.1");
            int i=0;
            for(pmptr pm2=pm1->next;pm2!=NULL;pm2=pm2->next){
                i++;
                //printf("\ninRep 2.11 %d",i);
                if(temp->utlz<pm2->utlz){
                    temp=pm2;
                }
            }
            if(temp==pm1){
                pm1=pm1->next;
                continue;
            }
            pmptr pmnext=pm1->next;

            switchpm(pm1,temp,t1);

            pm1=pmnext;
        }
    }
//printf("\ninRep 10");
}

int dbl_rep_heuristics(geneptr t1,pmptr pm1,int cn_ogi,int* cn_absent,int* absent_count){
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];

    float best_res=-1;
    cnptr cn11,cn22;
    vmptr vm11;
    for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                continue;
            }
            for(cnptr cn1=vm1->hcn; cn1!=NULL; cn1=cn1->next){
                float cpu_diff=vm1->cpu + cn1->cpu - cpu_req;
                float mem_diff=vm1->mem + cn1->mem - mem_req;

                if((cpu_diff<0)||(mem_diff<0)){
                    continue;
                }
                for(cnptr cn2=cn1->next; cn2!=NULL; cn2=cn2->next){

                    float cpu_rep=cn1->cpu+cn2->cpu;
                    float mem_rep=cn1->mem+cn2->mem;


                    if((cpu_rep+vm1->cpu-cpu_req<0)||(mem_rep+vm1->mem-mem_req<0)){
                        continue;
                    }
                    float vm_cpu=data1.VM_TYPE[vm1->ogi][0];
                    float vm_mem=data1.VM_TYPE[vm1->ogi][1];

                    float res_req=(cpu_req/vm_cpu)*(mem_req/vm_mem);
                    float res_rep=(cpu_rep/vm_cpu)*(mem_rep/vm_mem);
                    float res_diff=res_req-res_rep;

                    if(res_diff<=0){
                        continue;
                    }
                    if(res_diff>best_res){
                        vm11=vm1;
                        cn11=cn1;
                        cn22=cn2;
                        best_res=res_diff;
                    }
                }
            }


    }
    if(best_res>0){
        //printf("\n dbl rep done");
        //printf("\nAbsent Count:%d ",*absent_count);
        cn_absent[*absent_count-1]=cn11->ogi;
        cn_absent[*absent_count]=cn22->ogi;
        *absent_count=*absent_count+1;

        //printf("\nIN HEURISTICS 5");
        removeContainer(cn11->ogi,t1);
        removeContainer(cn22->ogi,t1);

        //printf("\nIN HEURISTICS 6");
        cnptr newCn=create_cn(cn_ogi);
        addCnInVm(newCn,vm11);
        //printf("\nIN HEURISTICS 7");
        //printf("\nAbsent Count2:%d ",*absent_count);
        return(1);
    }
    //printf("\n dbl rep not done");
    return(0);
}

int complete_fill_heuristics(pmptr pm1,int cn_ogi){
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];

    float cmplt_thrshld=.9;
    for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
        if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                continue;
        }
        // printf("\nIN HEURISTICS 3");
        float res_rem=(cpu_req/vm1->cpu)*(mem_req/vm1->mem);
        if(res_rem>cmplt_thrshld){
            //printf("\nIN HEURISTICS 4");
            cnptr cn1=create_cn(cn_ogi);
            addCnInVm(cn1,vm1);
            return(1);
        }
    }
    return(0);
}

int better_fill_heuristics(geneptr t1,pmptr pm1,int cn_ogi,int *cn_absent,int *absent_count){
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];

    float bttr_res=-1;
    vmptr vm11;
    cnptr cn11;
    for(vmptr vm1=pm1->hvm; vm1!=NULL; vm1=vm1->next){
            if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                continue;
            }
            for(cnptr cn1=vm1->hcn; cn1!=NULL; cn1=cn1->next){
                float cpu_diff=vm1->cpu + cn1->cpu - cpu_req;
                float mem_diff=vm1->mem + cn1->mem - mem_req;

                if((cpu_diff<0)||(mem_diff<0)){
                    continue;
                }
                float vm_cpu=data1.VM_TYPE[vm1->ogi][0];
                float vm_mem=data1.VM_TYPE[vm1->ogi][1];

                float res_req=(cpu_req/vm_cpu)*(mem_req/vm_mem);
                float res_rep=(cn1->cpu/vm_cpu)*(cn1->mem/vm_mem);
                float res_diff=res_req-res_rep;

                if(res_diff<=0){
                    continue;
                }
                if(res_diff>bttr_res){
                    bttr_res=res_diff;
                    cn11=cn1;
                    vm11=vm1;
                }

            }
    }

    if(bttr_res>0){
        //printf("\nAbsent Count2:%d ",*absent_count);
        cn_absent[*absent_count-1]=cn11->ogi;
        //printf("\nIN HEURISTICS 5");
        removeContainer(cn11->ogi,t1);
        //printf("\nIN HEURISTICS 6");
        cnptr newCn=create_cn(cn_ogi);
        addCnInVm(newCn,vm11);
        //printf("\nIN HEURISTICS 7");


        return(1);
    }

    return(0);
}

int best_fill_heuristics(pmptr pm1,int cn_ogi){
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];

    vmptr best_vm=NULL;
    float best_res=99999;
    //printf("\nIN HEURISTICS best fit");
    for(vmptr vm1=pm1->hvm; vm1!=NULL; vm1=vm1->next){
            if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                    continue;
            }
            float vm_cpu=data1.VM_TYPE[vm1->ogi][0];
            float vm_mem=data1.VM_TYPE[vm1->ogi][1];

            float res_rem=(((vm1->cpu-cpu_req)/vm_cpu)*((vm1->mem-mem_req)/vm_mem));
            if(res_rem<best_res){
                best_res=res_rem;
                best_vm=vm1;
            }
    }
    if(best_vm!=NULL){
        cnptr newCn=create_cn(cn_ogi);
        addCnInVm(newCn,best_vm);
        return(1);
    }
    return(0);
}

int first_fill_heuristics(pmptr pm1, int cn_ogi){
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];


    //printf("\nIN HEURISTICS best fit");
    for(vmptr vm1=pm1->hvm; vm1!=NULL; vm1=vm1->next){
            if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                    continue;
            }
            cnptr newCn=create_cn(cn_ogi);
            addCnInVm(newCn,vm1);
            return(1);
    }

    return(0);
}

void calculateContRes(geneptr t1){
    for(pmptr p1=t1->hpm;p1!=NULL;p1=p1->next){
        for(vmptr v1=p1->hvm;v1!=NULL;v1=v1->next){
            for(cnptr c1=v1->hcn;c1!=NULL;c1=c1->next){
                c1->norm_res=(c1->cpu/data1.pm[0])*(c1->mem/data1.pm[1]);
            }
        }
    }
}

cnptr findprevcn(cnptr c1,vmptr v1){
    for(cnptr cn1=v1->hcn;cn1->next!=NULL;cn1=cn1->next){
        if(cn1->next==c1){
            return(cn1);
        }
    }
    printf("ERROR:CANT FIND PREV CN");
    exit(0);
}

void swapContainers(cnptr c1,cnptr c2,vmptr v1){
    if((c1==NULL)||(c2==NULL)){
        printf("\nERROR: Cant swap null container.");
        exit(0);
    }
    if(c1==c2){
        return;
    }
    //if c1 and c2 are next to each other
    if(c1->next==c2){
        if(v1->hcn==c1){
            if(c2->next==NULL){
                v1->hcn=c2;
                c2->next=c1;
                c1->next=NULL;
                return;
            }
            cnptr c2next=c2->next;
            v1->hcn=c2;
            c2->next=c1;
            c1->next=c2next;
            return;
        }
        if(c2->next==NULL){
            cnptr c1prev=findprevcn(c1,v1);
            c1prev->next=c2;
            c2->next=c1;
            c1->next=NULL;
            return;
        }
        cnptr c1prev=findprevcn(c1,v1);
        cnptr c2next=c2->next;
        c1prev->next=c2;
        c2->next=c1;
        c1->next=c2next;
        return;
    }
    if(c2->next==c1){
        swapContainers(c2,c1,v1);
        return;
    }
    if(v1->hcn==c1){
        if(c2->next==NULL){
            cnptr c2prev=findprevcn(c2,v1);
            v1->hcn=c2;
            c2->next=c1->next;
            c2prev->next=c1;
            c1->next=NULL;
            return;
        }
        cnptr c2prev=findprevcn(c2,v1);
        cnptr c2next=c2->next;
        v1->hcn=c2;
        c2->next=c1->next;
        c2prev->next=c1;
        c1->next=c2next;
        return;
    }
    if(v1->hcn==c2){
        swapContainers(c2,c1,v1);
        return;
    }

    if(c2->next==NULL){
        cnptr c1prev=findprevcn(c1,v1);
        cnptr c2prev=findprevcn(c2,v1);
        c1prev->next=c2;
        c2->next=c1->next;
        c2prev->next=c1;
        c1->next=NULL;
        return;
    }
    if(c1->next==NULL){
        swapContainers(c2,c1,v1);
        return;
    }

    cnptr c1prev=findprevcn(c1,v1);
    cnptr c2prev=findprevcn(c2,v1);
    cnptr c2next=c2->next;
    c1prev->next=c2;
    c2->next=c1->next;
    c2prev->next=c1;
    c1->next=c2next;
    return;
}

void sortContainersByRes(geneptr t1){
    calculateContRes(t1);
    for(pmptr p1=t1->hpm;p1!=NULL;p1=p1->next){
        for(vmptr v1=p1->hvm;v1!=NULL;v1=v1->next){
            for(cnptr c1=v1->hcn;c1!=NULL;c1=c1->next){
                cnptr best_cn=c1;
                for(cnptr c2=c1->next;c2!=NULL;c2=c2->next){
                    if(best_cn->norm_res>c2->norm_res){
                        best_cn=c2;
                    }
                }
                if(best_cn==c1){
                    continue;
                }
                cnptr tempcn=best_cn;
                swapContainers(c1,best_cn,v1);
                c1=tempcn;
            }
        }
    }
}

void ff_rc_heuristics(cnptr cn,geneptr t1){
    float cpu_req=cn->cpu;
    float mem_req=cn->mem;
    cnptr new_cn=create_cn(cn->ogi);

    vmptr selected_vm=NULL;
    for(pmptr p1=t1->hpm;p1!=NULL;p1=p1->next){
        for(vmptr v1=p1->hvm;v1!=NULL;v1=v1->next){
            if((v1->cpu>=cpu_req)&&(v1->mem>=mem_req)){
                selected_vm=v1;
                break;
            }
        }
    }
    //if theres any vm that can contain both
    if(selected_vm!=NULL){
        addCnInVm(new_cn,selected_vm);
        return;
    }

    //create new vm
    int vm_ogi=rand()%VM_SIZE;
    float vm_cpu=data1.VM_TYPE[vm_ogi][0];
    float vm_mem=data1.VM_TYPE[vm_ogi][1];

    while((data1.pm[0]<(vm_cpu*1.1)) || (data1.pm[1]<(vm_mem+200)) || (vm_mem<mem_req) || (vm_cpu<cpu_req)){
        vm_ogi=rand()%VM_SIZE;
        vm_cpu=data1.VM_TYPE[vm_ogi][0];
        vm_mem=data1.VM_TYPE[vm_ogi][1];
    }

    vmptr new_vm=create_vm(vm_ogi);
    addCnInVm(new_cn,new_vm);

    //find pm which has space for this vm
    pmptr selected_pm=NULL;
    float req_cpu_vm=vm_cpu*1.1;
    float req_mem_vm=vm_mem+200;

    for(pmptr p1=t1->hpm;p1!=NULL;p1=p1->next){
        if((p1->cpu>=req_cpu_vm)&&(p1->mem>=req_mem_vm)){
            selected_pm=p1;
            break;
        }
    }
    //if there is an existing pm that can have it
    if(selected_pm!=NULL){
        addVmInPm(new_vm,selected_pm);
        return;
    }
    //if there's not
    pmptr new_pm=create_pm(t1);
    addVmInPm(new_vm,new_pm);

}

int rep_small_two_heuristics(geneptr t1,vmptr v1,int cn_absent[CNT_SIZE],int absent_count){

    for(int i=0;i<absent_count;i++){
        float cpu_req=data1.CNT_TYPE[cn_absent[i]][0];
        float mem_req=data1.CNT_TYPE[cn_absent[i]][1];


        cnptr first_cn=v1->hcn;
        cnptr second_cn=first_cn->next;

        if((first_cn->cpu+second_cn->cpu<=cpu_req)&&(first_cn->mem+second_cn->mem<=mem_req)){
            //replace the two cn with given cn
                //delete 2 cn
            removeContainer(first_cn->ogi,t1);
            removeContainer(second_cn->ogi,t1);
                //add new container to vm
            cnptr new_cn=create_cn(cn_absent[i]);
            addCnInVm(new_cn,v1);

            //allocate the two cn with ff/ff&RC. create a random vm if no space is available in any vm.
            ff_rc_heuristics(first_cn,t1);
            ff_rc_heuristics(second_cn,t1);

            //sort all again
            sortContainersByRes(t1);
            //return the id of 2 containers deleted.
            return(i);
        }
    }
    return(-1);
}

void allocate_list_ff_rc(geneptr t1,int cn_absent[CNT_SIZE],int absent_count){
    //loop over all one by one
    for(int i=0;i<absent_count;i++){
        //allocate using ff_rc/ff
        cnptr newcn=create_cn(cn_absent[i]);
        ff_rc_heuristics(newcn,t1);
    }

}

void heuristics2(geneptr t1,int cn_absent[CNT_SIZE],int absent_count){
    //if(absent_count==0){
     //   return;
    //}
    //sortPmByUtlz(t1,1);
    //printf("\nAbsent Count:%d ",absent_count);
    //int cn_ogi=cn_absent[absent_count-1];
    //float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    //float mem_req=data1.CNT_TYPE[cn_ogi][1];

    sortContainersByRes(t1);
    //printstats(t1);

    for(pmptr p1=t1->hpm;p1!=NULL;p1=p1->next){
        vmptr v1=p1->hvm;
        while((v1!=NULL)&&(absent_count>0)){
            int replaced_container=rep_small_two_heuristics(t1,v1,cn_absent,absent_count);

            if(replaced_container==-1){
                v1=v1->next;
                continue;
            }
            if(absent_count==(replaced_container+1)){
                absent_count--;
                continue;
            }
            cn_absent[replaced_container]=cn_absent[absent_count-1];
            absent_count--;
        }
    }

    if(absent_count>0){
        allocate_list_ff_rc(t1,cn_absent,absent_count);
    }

}

void heuristics(geneptr t1,int cn_absent[CNT_SIZE],int absent_count){
    if(absent_count==0){
        return;
    }
    sortPmByUtlz(t1,1);
    //printf("\nAbsent Count:%d ",absent_count);
    int cn_ogi=cn_absent[absent_count-1];
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];


    //double replacement.
     //cnptr cn1db,cn2db;
    //vmptr vm1db;


    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){

        int dc1=dbl_rep_heuristics(t1,pm1,cn_ogi,cn_absent,&absent_count);
        if(dc1==1){
            //printf("\nAbsent Count here:%d ",absent_count);
            heuristics(t1,cn_absent,absent_count);
            return;

        }
        int cmp_fill=complete_fill_heuristics(pm1,cn_ogi);
        if(cmp_fill){
            absent_count--;
            heuristics(t1,cn_absent,absent_count);
            return;
        }
        int bttr_fill=better_fill_heuristics(t1,pm1,cn_ogi,cn_absent,&absent_count);
        if(bttr_fill){
            heuristics(t1,cn_absent,absent_count);
            return;
        }
        int ff=first_fill_heuristics(pm1,cn_ogi);
        if(ff){
            absent_count--;
            heuristics(t1,cn_absent,absent_count);
            return;
        }
    }

    //cant allocate to any existing vm/pm, create a new one

    int vm_ogi=ff_vm_new(cn_ogi);

    vmptr newVm=create_vm(vm_ogi);
    cnptr newCn=create_cn(cn_ogi);

    pmptr pm2;
    int foundPm=0;
    //finding any existing pms that can contain new vm.
    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        int hasSpace=(pm1->cpu>(newVm->cpu*1.1)) && (pm1->mem>(newVm->mem+200));

        if(hasSpace){
            foundPm=1;
            pm2=pm1;
            break;
        }
    }
    if(!foundPm){
        //printstats(t1);
        pm2=create_pm(t1);
    }
    addVmInPm(newVm,pm2);
    addCnInVm(newCn,newVm);
    fitness1(t1);

    absent_count--;
    heuristics(t1,cn_absent,absent_count);
    return;

}

geneptr makeReplicaGene(geneptr t1){
    geneptr t2=create_gene();
    //printf("\ninRep 1");
    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        pmptr pm2=create_pm(t2);
        for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            vmptr vm2=create_vm(vm1->ogi);
            for(cnptr cn1=vm1->hcn; cn1!=NULL; cn1=cn1->next){
                    cnptr cn2=create_cn(cn1->ogi);
                    addCnInVm(cn2,vm2);
            }
            addVmInPm(vm2,pm2);
        }
    }
    fitness1(t2);
    //printf("\ninRep 10");
    return(t2);
}

geneptr mutation(popptr p1){
    geneptr parent1=binary_trnment(p1);
    fitness1(parent1);
    parent1=makeReplicaGene(parent1);

    int remove_n=MUTATION_VM_REM;
    int absent_cn[CNT_SIZE];
    int absent_count=0;
    for(pmptr pm1=parent1->hpm;pm1!=NULL;pm1=pm1->next){
        if(remove_n==0){
            break;
        }
        for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            if(vm1->utlz>MUTATION_VM_THRSHLD){
                continue;
            }
            int dc=rand()%10;
            if(dc>=MUTATION_VM_PROB){
                continue;
            }
            for(cnptr cn1=vm1->hcn; cn1!=NULL; cn1=cn1->next){
                absent_cn[absent_count]=cn1->ogi;
                absent_count++;
            }
            remove_vm(pm1,vm1);
            remove_n--;
        }

    }

    heuristics(parent1,absent_cn,absent_count);
    return(parent1);

}

void crossMut(popptr p1){
    int prob_crs=rand()%10;
    geneptr child;
    if(prob_crs<PROB_CROSSOVER){
        child=crossover(p1);
    }
    else
        child=mutation(p1);

    fitness1(child);
    unpack_pm(child);
    fitness(child);
    merge_vm(child);
    printstats(child);
    fitness(child);
    add_gene(p1,child);
}

geneptr crossover2(geneptr parent1,geneptr parent2,int cpumem){
    geneptr child=create_gene();
    //int min_pm=parent1->npm>parent2->npm?parent2->npm:parent1->npm;
    parent1=makeReplicaGene(parent1);
    parent2=makeReplicaGene(parent2);

    sortPmByUtlz(parent1,cpumem);
    sortPmByUtlz(parent2,cpumem);
    int i=0;
    pmptr parent1_cpm=parent1->hpm;
    pmptr parent2_cpm=parent2->hpm;

    while((parent1_cpm!=NULL)||(parent2_cpm!=NULL)){
        pmptr best_pm;
        geneptr best_parent,othr_parent;
        if(cpumem==0){
            best_pm=parent1_cpm->mem_utlz>parent2_cpm->mem_utlz?parent1_cpm:parent2_cpm;
            best_parent=parent1_cpm->mem_utlz>parent2_cpm->mem_utlz?parent1:parent2;
            othr_parent=best_parent==parent1?parent2:parent1;
        }
        else{
            best_pm=parent1_cpm->utlz>parent2_cpm->utlz?parent1_cpm:parent2_cpm;
            best_parent=parent1_cpm->utlz>parent2_cpm->utlz?parent1:parent2;
            othr_parent=best_parent==parent1?parent2:parent1;
        }

        pmptr cpm=create_pm(child);
        copy_pm(best_pm,cpm);
        printf("\nin crsover 4");
        //remove pm's containers in other parent.

        vmptr vm1=best_pm->hvm;
        while(vm1!=NULL){
            cnptr cn1=vm1->hcn;
            while(cn1!=NULL){
                removeContainer(cn1->ogi,othr_parent);
                cn1=cn1->next;
            }
            vm1=vm1->next;
        }


        remove_pm(best_pm,best_parent);

        parent1_cpm=parent1_cpm->next;
        parent2_cpm=parent2_cpm->next;
    }
    printf("\nin crsover 5");
    int cn_present[CNT_SIZE];
    count_containers(child);
    cnInGene(child,cn_present);

    int cn_absent[CNT_SIZE];
    int absent_count=0;
    for(int i=0;i<CNT_SIZE;i++){
        if(cn_present[i]==0){
            cn_absent[absent_count]=i;
            absent_count++;
        }
    }

    //check the integrity of the list.1.no duplicates checked


    if(absent_count!=0){
        printf("in crossover2 7");
        heuristics2(child,cn_absent,absent_count);
    }

    count_containers(child);
    integrity_duplicate_containers(child);
    integrity_all_containers(child);

    return child;
}

geneptr crossover(popptr p1){
    geneptr parent1=binary_trnment(p1);
    geneptr parent2=binary_trnment(p1);

    while(parent1==parent2){
        parent2=binary_trnment(p1);
    }

    printf("\n\nparent 1 utilization:%f,parent 2 utilization:%f",parent1->utlz,parent2->utlz);

    geneptr child1=crossover2(parent1, parent2, 0);
    geneptr child2=crossover2(parent1, parent2, 1);

    return(child1);
}

pmptr search_best_pm(geneptr t1){
    if(t1->hpm==NULL){
        return(NULL);
    }
    pmptr pm=t1->hpm,bestpm;
    float best_utlz=0;
    while(pm!=NULL){
        if(best_utlz<pm->utlz){
            best_utlz=pm->utlz;
            bestpm=pm;
        }
        pm=pm->next;
    }
    return bestpm;
}

geneptr binary_trnment(popptr p1){
    int i1=rand()%p1->n;
    int i2=rand()%p1->n;

    while(i2==i1){
        i2=rand()%p1->n;
    }

    geneptr temp,t1,t2;

    t1=p1->gene;
    t2=p1->gene;

    for(int i=0;i<i1;i++){
        t1=t1->next;
    }
    for(int i=0;i<i2;i++){
        t2=t2->next;
    }
    if(t1->utlz<t2->utlz){
        temp=t1;
        t1=t2;
        t2=temp;
    }

    int dc=rand()%10;
    if(dc<PROB_GENE_BIN_TRN){
        return t1;
    }
    return t2;
}

void readCsv(){

    FILE *fpr;
    FILE *fpr2;
    FILE *fpr3;
    FILE *fpr4;
    FILE *fpr5;
    FILE *fpr6;
    FILE *fpr7;
    FILE *fpr8;
    FILE *fpr9;
    FILE *fpr10;
    FILE *fpr11;

    char c;
    fpr=fopen("data/Container200_ten.csv","r");
    fpr2=fopen("data/Container200_twenty.csv","r");
    fpr3=fopen("data/VMConfig_ten.csv","r");
    fpr4=fopen("data/VMConfig_twenty.csv","r");
    fpr5=fopen("data/PMConfig_small.csv","r");
    fpr6=fopen("data/Container500_ten.csv","r");
    fpr7=fopen("data/Container500_twenty.csv","r");
    fpr8=fopen("data/Container1000_ten.csv","r");
    fpr9=fopen("data/Container1000_twenty.csv","r");
    fpr10=fopen("data/Container1500_ten.csv","r");
    fpr11=fopen("data/Container1500_twenty.csv","r");


    if(fpr==NULL||fpr2==NULL||fpr3==NULL||fpr4==NULL||fpr5==NULL||fpr6==NULL||fpr7==NULL||fpr8==NULL||fpr9==NULL||fpr10==NULL||fpr11==NULL){
        printf("\n ERROR in opening file\n");
    }
    //reading container file for synthetic vms
    int i=0;
    while(i<200){
        fscanf(fpr,"%f,%f\n",&data1.cont_syn200[i][0],&data1.cont_syn200[i][1]);
        i++;
    }
    //reading container file for real vms
    i=0;
    while(i<200){
        fscanf(fpr2,"%f,%f\n",&data1.cont_rea200[i][0],&data1.cont_rea200[i][1]);
        i++;
    }
    i=0;

    while(i<200){
        fscanf(fpr6,"%f,%f\n",&data1.cont_syn500[i][0],&data1.cont_syn500[i][1]);
        i++;
    }
    //reading container file for real vms
    i=0;
    while(i<200){
        fscanf(fpr7,"%f,%f\n",&data1.cont_rea500[i][0],&data1.cont_rea500[i][1]);
        i++;
    }

    i=0;

    while(i<200){
        fscanf(fpr8,"%f,%f\n",&data1.cont_syn1000[i][0],&data1.cont_syn1000[i][1]);
        i++;
    }
    //reading container file for real vms
    i=0;
    while(i<200){
        fscanf(fpr9,"%f,%f\n",&data1.cont_rea1000[i][0],&data1.cont_rea1000[i][1]);
        i++;
    }

    i=0;

    while(i<200){
        fscanf(fpr10,"%f,%f\n",&data1.cont_syn1500[i][0],&data1.cont_syn1500[i][1]);
        i++;
    }
    //reading container file for real vms
    i=0;
    while(i<200){
        fscanf(fpr11,"%f,%f\n",&data1.cont_rea1500[i][0],&data1.cont_rea1500[i][1]);
        i++;
    }


    //reading synthetic vm config
    i=0;
    while(i<10){
        fscanf(fpr3,"%f,%f\n",&data1.vm_syn[i][0],&data1.vm_syn[i][1]);
        i++;
    }

    //reading real vm config
    i=0;
    while(i<20){
        fscanf(fpr4,"%f,%f\n",&data1.vm_rea[i][0],&data1.vm_rea[i][1]);
        i++;
    }

    //reading pm config
    fscanf(fpr5,"%f\n%f\n%f",&data1.pm[0],&data1.pm[1],&data1.pm[2]);
    fclose(fpr);
    fclose(fpr2);
    fclose(fpr3);
    fclose(fpr4);
    fclose(fpr5);

}

void fitness1(geneptr t1){
    t1->utlz=0;
    pmptr pm1=t1->hpm;
    int i=0;
    //printf("\n\nCalculating fitness...");
    while(pm1!=NULL){
        vmptr vm1=pm1->hvm;
        float utlz=0;
        float mem_u=0;
        while(vm1!=NULL){
            float vm_utlz=0;
            float vm_mem_utlz=0;
            cnptr cn1=vm1->hcn;
            while(cn1!=NULL){
                utlz+=cn1->cpu;
                mem_u+=cn1->mem;
                vm_utlz+=cn1->cpu;
                vm_mem_utlz+=cn1->mem;
                cn1=cn1->next;
            }
            int ogi=vm1->ogi;
            float vm_overhead=data1.VM_TYPE[ogi][0]*0.1;
            float vm_mem=data1.VM_TYPE[ogi][1];
            vm_utlz=vm_utlz/(vm_overhead*10);
            vm_mem_utlz=vm_mem_utlz/vm_utlz;
            mem_u+=200;
            vm1->utlz=vm_utlz;
            utlz+=vm_overhead;
            vm1=vm1->next;
        }
        utlz=utlz/data1.pm[0];
        pm1->utlz=utlz;
        pm1->mem_utlz=mem_u/data1.pm[1];
        t1->utlz+=utlz;
       // printf("\n\tPM:%dCpu Utlz:%f,\tMem Utlz:%f",i,pm1->utlz,pm1->mem_utlz);
        pm1=pm1->next;
        i++;
    }
    t1->utlz=t1->utlz/t1->npm;
}

void fitness(geneptr t1){
    t1->utlz=0;
    pmptr pm1=t1->hpm;
    int i=0;
    printf("\n\nCalculating fitness...");
    while(pm1!=NULL){
        vmptr vm1=pm1->hvm;
        float utlz=0;
        float mem_u=0;
        while(vm1!=NULL){
            float vm_utlz=0;
            float vm_mem_utlz=0;
            cnptr cn1=vm1->hcn;
            while(cn1!=NULL){
                utlz+=cn1->cpu;
                mem_u+=cn1->mem;
                vm_utlz+=cn1->cpu;
                vm_mem_utlz+=cn1->mem;
                cn1=cn1->next;
            }
            int ogi=vm1->ogi;
            float vm_overhead=data1.VM_TYPE[ogi][0]*0.1;
            float vm_mem=data1.VM_TYPE[ogi][1];
            vm_utlz=vm_utlz/(vm_overhead*10);
            vm_mem_utlz=vm_mem_utlz/vm_utlz;
            mem_u+=200;
            vm1->utlz=vm_utlz;
            utlz+=vm_overhead;
            vm1=vm1->next;
        }
        utlz=utlz/data1.pm[0];
        pm1->utlz=utlz;
        pm1->mem_utlz=mem_u/data1.pm[1];
        t1->utlz+=utlz;
        printf("\n\tPM:%dCpu Utlz:%f,\tMem Utlz:%f",i,pm1->utlz,pm1->mem_utlz);
        pm1=pm1->next;
        i++;
    }
    t1->utlz=t1->utlz/t1->npm;

}

void createPop(){
    head_population=(popptr)malloc(sizeof(struct population));
    head_population->n=0;
    popptr h1=head_population;

    //geneptr t1;
    //geneptr t2=(geneptr)malloc(sizeof(struct gene));
    //h1->gene=t2;
    //printf("\ncreatepop 1");
    //t2->npm=0;
    //t2->hpm=NULL;
    for(int i=0;i<NO_POP;i++){
        //t1=t2;
        geneptr t1=create_gene();
        //printf("\ncreatepop 2");
        generate_random_containers(t1);
        //printf("\ncreatepop 3");
        first_fit(t1);
       // printf("\ncreatepop 4");
       // printstats(t1);
       // printf("\ncreatepop 5");
        fitness(t1);
        integrity_all_containers(t1);
        add_gene(h1,t1);
        //t2=(geneptr)malloc(sizeof(struct gene));
        //t1->next=t2;
        //t2->npm=0;
        //t2->hpm=NULL;
    }

    search_best_gene(h1);
}

void search_best_gene(popptr p1){
    geneptr t1=p1->gene;
    float best_utlz=0;
    int i=0;
    printf("\n\nFinding gene with best Utlz");
    while(t1!=NULL){
        printf("\n\tGene:%d Utlz:%f",i,t1->utlz);
        if(t1->utlz>best_utlz){
            best_utlz=t1->utlz;
        }

        t1=t1->next;
        i++;
    }
    printf("\nBest Utlz:%f",best_utlz);
}

void generate_random_containers(geneptr t1){
    int i=0;
    int n=CNT_SIZE;
    int check_bucket[n];
    while(i<n){
        check_bucket[i]=0;
        i++;
    }
    i=0;
    while(i<n){
        int ran=(rand()%n);
        while(check_bucket[ran]!=0){
            ran=(rand()%n);
        }
        t1->test_containers[i]=ran;
        check_bucket[ran]=1;
        //printf("\n%d added in test_containers\n",ran);
        i++;
    }

}

void printstats(geneptr t1){
    pmptr pm=t1->hpm;
    int i=0;
    while(pm!=NULL){
        printf("\n\nPM:%d : CPU :%f/%f Mem:%f/%f No of VM's:%d",i,pm->cpu,data1.pm[0],pm->mem,data1.pm[1],pm->nvm);
        int j=0;
        vmptr vm=pm->hvm;
        while(vm!=NULL){
            printf("\n\tVM:%d : CPU:%f/%f Mem:%f/%f No of Containers:%d Utlz:%f",j,vm->cpu,data1.VM_TYPE[vm->ogi][0],vm->mem,data1.VM_TYPE[vm->ogi][1],vm->ncn,vm->utlz);
            cnptr cn=vm->hcn;
            int k=0;
            while(cn!=NULL){
                printf("\n\t\tContainer:%d : %d : CPU taken:%f Mem taken:%f",k,cn->ogi,cn->cpu,cn->mem);
                cn=cn->next;
                k++;
            }

            vm=vm->next;
            j++;
        }
        pm=pm->next;
        i++;
    }

}


//--first fit--
int bf_vm_new(int cn_ogi){
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];

    float best_res=-1;
    int vmi=-1;
    for(int i=0;i<VM_SIZE;i++){
        int vm_ogi=i;
        float vm_cpu=data1.VM_TYPE[vm_ogi][0];
        float vm_mem=data1.VM_TYPE[vm_ogi][1];

        if((data1.pm[0]<(vm_cpu*1.1)) || (data1.pm[1]<(vm_mem+200)) || (vm_mem<mem_req) || (vm_cpu<cpu_req)){
            continue;
        }
        float res_diff=(vm_cpu-cpu_req)-(vm_mem-mem_req);
        res_diff=abs(res_diff);
        if(best_res<res_diff){
            best_res=res_diff;
            vmi=i;
        }

    }
    if(vmi==-1){
        printf("\nERROR:NO VM CAN CONTAIN THIS CONTAINER");
        exit(0);
    }
    return(vmi);
}

int ff_vm_new(int cn_ogi){
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];

    int vm_ogi=rand()%VM_SIZE;
    float vm_cpu=data1.VM_TYPE[vm_ogi][0];
    float vm_mem=data1.VM_TYPE[vm_ogi][1];

    while((data1.pm[0]<(vm_cpu*1.1)) || (data1.pm[1]<(vm_mem+200)) || (vm_mem<mem_req) || (vm_cpu<cpu_req)){
        vm_ogi=rand()%VM_SIZE;
        vm_cpu=data1.VM_TYPE[vm_ogi][0];
        vm_mem=data1.VM_TYPE[vm_ogi][1];
    }
    return(vm_ogi);
}

vmptr ff_vm(geneptr t1,float cpureq,float memreq){
        for(int j=0;j<t1->npm;j++){
            pmptr pm=t1->hpm;
            while(pm!=NULL){
                vmptr vm=pm->hvm;
                pm=pm->next;
                while(vm!=NULL){
                    float vmcpu=vm->cpu;
                    float vmmem=vm->mem;

                    if((cpureq>vmcpu)||(memreq>vmmem)){
                        vm=vm->next;
                        continue;
                    }
                    return(vm);
                }
            }
        }
        return NULL;
}
pmptr ff_pm(geneptr t1,float vmCpu,float vmMem){
for(int j=0;j<t1->npm;j++){
    pmptr pm=t1->hpm;
    while(pm!=NULL){
        float pmCpu=pm->cpu;
        float pmMem=pm->mem;


        if((pmCpu<(vmCpu*1.1))||(pmMem<(vmMem+200))){
            pm=pm->next;
            continue;
        }
        return(pm);


        }
    }
    return(NULL);

}
pmptr create_pm(geneptr t1){
    if(t1->npm==0){
        pmptr pm=(pmptr)malloc(sizeof(struct pm));
        pm->cpu=data1.pm[0];
        pm->mem=data1.pm[1];
        pm->next=NULL;
        pm->utlz=0;
        pm->nvm=0;
        pm->hvm=NULL;
        t1->hpm=pm;
        t1->npm++;
        return(pm);
    }
    pmptr pm=t1->hpm;
    while(pm->next!=NULL){
        pm=pm->next;
    }
    pmptr pm1=(pmptr)malloc(sizeof(struct pm));
    pm1->cpu=data1.pm[0];
    pm1->mem=data1.pm[1];
    pm1->next=NULL;
    pm1->utlz=0;
    pm1->hvm=NULL;
    pm->next=pm1;
    t1->npm++;
    return(pm1);

}

vmptr create_vm(int vmi){
    vmptr vm=(vmptr)malloc(sizeof(struct vm));
    vm->cpu=data1.VM_TYPE[vmi][0];
    vm->mem=data1.VM_TYPE[vmi][1];

    vm->hcn=NULL;
    vm->ncn=0;
    vm->next=NULL;
    vm->ogi=vmi;
    vm->utlz=0;

    return(vm);

}
cnptr create_cn(int cni){
    cnptr cn=(cnptr)malloc(sizeof(struct cn));
    cn->cpu=data1.CNT_TYPE[cni][0];
    cn->mem=data1.CNT_TYPE[cni][1];
    cn->next=NULL;
    cn->ogi=cni;
    return(cn);
}
void addVmInPm(vmptr vm,pmptr pm){
    float vm_cpu=data1.VM_TYPE[vm->ogi][0];
    float vm_mem=data1.VM_TYPE[vm->ogi][1];

    pm->cpu-=(vm_cpu*1.1);
    pm->mem-=(vm_mem+200);
    if((pm->cpu<0)||(pm->mem<0)){
        printf("\nPm space integrity error\n");
        exit(0);
    }

    if(pm->hvm==NULL){
        pm->hvm=vm;
    }
    else{
        //printf("\naddvminpm 4");
        vmptr vm2=pm->hvm;

        while(vm2->next!=NULL){
            vm2=vm2->next;
        }

        //printf("\naddvminpm 6");
        vm2->next=vm;
    }
    //printf("\naddvminpm 7");
    pm->nvm+=1;
}

void addCnInVm(cnptr cn,vmptr vm){
    //printf("\naddcn:mem:%f",vm->mem);
    vm->cpu-=cn->cpu;
    vm->mem-=cn->mem;
   // printf("\naddcn:mem:%f",vm->mem);
    if((vm->cpu<0)||(vm->mem<0)){
        printf("\nVm space integrity error,%f - %f,%f - %f\n",vm->cpu+cn->cpu,cn->cpu,vm->mem+cn->cpu,cn->mem);
        exit(0);
    }
    if(vm->ncn==0){
        vm->hcn=cn;
    }
    else{
        cnptr cn1=vm->hcn;
        while(cn1->next!=NULL){
            cn1=cn1->next;
        }
        cn1->next=cn;
    }
    vm->ncn++;
}

void first_fit(geneptr t1){
    printf("\nff 1");
    for(int i=0;i<CNT_SIZE;i++){
        //printf("\nff 2");
        float cpureq=data1.CNT_TYPE[t1->test_containers[i]][0];
        float memreq=data1.CNT_TYPE[t1->test_containers[i]][1];

        cnptr selected_cn=create_cn(t1->test_containers[i]);
        //printf("\nff 3");
        //if theres space in any vm
        vmptr selected_vm=ff_vm(t1,cpureq,memreq);
        //printf("\nff 4");



        if(selected_vm!=NULL){
            addCnInVm(selected_cn,selected_vm);
            //printf("\nff 5");
            continue;
        }

        int rndVm=ff_vm_new(t1->test_containers[i]);
        float vmCpu=data1.VM_TYPE[rndVm][0];
        float vmMem=data1.VM_TYPE[rndVm][1];

        pmptr selected_pm=ff_pm(t1,vmCpu,vmMem);

        if(selected_pm==NULL){
            selected_pm=create_pm(t1);
        }

        selected_vm=create_vm(rndVm);
        addVmInPm(selected_vm,selected_pm);
        //printf("\nff 9:%f -%f",vmMem,memreq);
        fflush(stdout);
        addCnInVm(selected_cn,selected_vm);
        //printf("\nff 10");
    }
}
