#include<stdio.h>
#include<stdlib.h>
#include<time.h>


#define NO_POP 10
#define NO_CROSSOVER 5
#define CNT_TYPE cont_syn200
#define CNT_SIZE 200
#define VM_TYPE vm_syn
#define VM_SIZE 10
#define UNPACK_THRSHLD .2
#define PROB_PM_TO_ADD 8
#define PROB_GENE_BIN_TRN 8

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

};
struct vm{
    int ogi;
    int ncn;
    float cpu;
    float mem;
    float utlz;
    struct vm *next;
    struct cn *hcn;
};
struct pm{
    int nvm;
    float cpu;
    float mem;
    float utlz;
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



void readCsv();
void createPop();
void generate_random_containers(geneptr t1);
void first_fit(geneptr t1);
void printstats(geneptr t1);
vmptr ff_vm(geneptr t1,float cpureq,float memreq);
pmptr ff_pm(geneptr t1,float vmCpu,float vmMem);
pmptr create_pm(geneptr t1);
vmptr create_vm(int vmi);
cnptr create_cn(int cni);
void addVmInPm(vmptr vm,pmptr pm);
void addCnInVm(cnptr cn,vmptr vm);
void printstats(geneptr t1);
void fitness(geneptr t1);
void search_best_gene(popptr p1);
void crossover(popptr p1);
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
void unpack_pm(geneptr t1,int* cn_present);
void merge_vm(geneptr t1);
void remove_vm(pmptr pm,vmptr vm);
void merge2vm(vmptr vm1,vmptr vm2,pmptr pm);

int main (){

    srand(time(0));
    readCsv();
    printf("\nread done");
    createPop();
    //printf("\ntotal %d genes",head_population->total);
    //for(int i=0;i<NO_CROSSOVER;i++){
   //     crossover();
    //}
   // search_best_gene();
    //printf("\ntotal %d genes",head_population->total);
    printf("\nend\n");
    for(int i=0;i<NO_CROSSOVER;i++){
        crossover(head_population);
    }
    search_best_gene(head_population);
}

void unpack_pm(geneptr t1,int* cn_present){
    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        if(pm1->utlz>UNPACK_THRSHLD){
            continue;
        }
        for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            for(cnptr cn1=vm1->hcn;cn1!=NULL;cn1=cn1->next){
                cn_present[cn1->ogi]=0;
            }
        }
        remove_pm(pm1,t1);
    }
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

void merge2vm(vmptr vm1,vmptr vm2,pmptr pm){
    //merge vm1 into vm2
    for(cnptr cn1=vm1->hcn;cn1!=NULL;cn1=cn1->next){
        cnptr newCn=create_cn(cn1->ogi);
        addCnInVm(newCn,vm2);
    }
    remove_vm(pm,vm1);
}

void remove_vm(pmptr pm,vmptr vm){
    //empty case
    if(pm->nvm==0){
        printf("\nERROR:NO VM FOUND TO REMOVE");
        exit(0);
    }
    //first value case
    vmptr vm2=pm->hvm;
    float vm_cpu=data1.VM_TYPE[vm->ogi][0];
    float vm_mem=data1.VM_TYPE[vm->ogi][1];
    pm->cpu+=(vm_cpu*1.1);
    pm->mem+=(vm_mem+200);


    if(vm2==vm){
        pm->hvm=vm2->next;
        pm->nvm--;
        return;
    }
    //any value
    while(vm2->next!=vm){
        vm2=vm2->next;
    }
    vm2->next=vm->next;
    pm->nvm--;
}

void merge_vm(geneptr t1){
    int count=0;
    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        int dc=1;
        vmptr vm1;
        //empty case
        if(pm1->hvm==NULL){
            continue;
        }
        while(dc){
            dc=0;
            vm1=pm1->hvm;
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
                dc=1;
                count++;
                break;
            }
        }
    }
    printf("\n%d vm merged",count);
    fitness(t1);
}

void heuristics(geneptr t1,int cn_absent[CNT_SIZE],int absent_count){
    if(absent_count==0){
        return;
    }
    //printf("\nAbsent Count:%d ",absent_count);
    int cn_ogi=cn_absent[absent_count-1];
    float cpu_req=data1.CNT_TYPE[cn_ogi][0];
    float mem_req=data1.CNT_TYPE[cn_ogi][1];

    pmptr best_pm;
    float best_utlz=-1;

    for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
        int hasSpace=0;
        for(vmptr vm1=pm1->hvm;vm1!=NULL;vm1=vm1->next){
            if((vm1->cpu>cpu_req)&&(vm1->mem>mem_req)){
                hasSpace=1;
                break;
            }
        }
        if(hasSpace){
            if(best_utlz<pm1->utlz){
                best_utlz=pm1->utlz;
                best_pm=pm1;
            }
        }
    }
    if(best_utlz==-1){

        int vm_ogi=bf_vm_new(cn_ogi);

        vmptr newVm=create_vm(vm_ogi);
        cnptr newCn=create_cn(cn_ogi);

        //printf("\nNO SPACE :%f\t,%f",newVm->cpu*1.1,mem_req,newVm->mem+200);
        pmptr pm2;
        int foundPm=0;
        for(pmptr pm1=t1->hpm;pm1!=NULL;pm1=pm1->next){
            int hasSpace=(pm1->cpu>(newVm->cpu*1.1)) && (pm1->mem>(newVm->mem+200));

            if(hasSpace){
                foundPm=1;
                //printf("\n\n\n\n\nhere");
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
        fitness(t1);

        absent_count--;
        heuristics(t1,cn_absent,absent_count);
        return;
    }
    //printf("\nIN HEURISTICS 2");
    //fill flags
    //int cmplt=0,bttr=0,best;
    //complete fill
    float cmplt_thrshld=10;
    for(vmptr vm1=best_pm->hvm;vm1!=NULL;vm1=vm1->next){
        if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                continue;
        }
       // printf("\nIN HEURISTICS 3");
        float cpu_rem=vm1->cpu-cpu_req;
        if(cpu_rem<cmplt_thrshld){
            //printf("\nIN HEURISTICS 4");
            cnptr cn1=create_cn(cn_ogi);
            addCnInVm(cn1,vm1);

            absent_count--;
            heuristics(t1,cn_absent,absent_count);
            return;
        }
    }

    //bttr fill
    for(vmptr vm1=best_pm->hvm; vm1!=NULL; vm1=vm1->next){
        if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                continue;
        }
        for(cnptr cn1=vm1->hcn; cn1!=NULL; cn1=cn1->next){
            if(cn1->cpu>=cpu_req){
                continue;
            }
            float cpu_diff=vm1->cpu + cn1->cpu - cpu_req;
            float mem_diff=vm1->mem + cn1->mem - mem_req;

            if((cpu_diff<0)||(mem_diff<0)){
                continue;
            }
            cn_absent[absent_count-1]=cn1->ogi;
            //printf("\nIN HEURISTICS 5");
            removeContainer(cn1->ogi,t1);
            //printf("\nIN HEURISTICS 6");
            cnptr newCn=create_cn(cn_ogi);
            addCnInVm(newCn,vm1);
            //printf("\nIN HEURISTICS 7");

            heuristics(t1,cn_absent,absent_count);
            return;
        }
    }

    //best_fit
    vmptr best_vm;
    float best_cpu=99999;
    //printf("\nIN HEURISTICS best fit");
    for(vmptr vm1=best_pm->hvm; vm1!=NULL; vm1=vm1->next){
        if((vm1->cpu<cpu_req)||(vm1->mem<mem_req)){
                continue;
        }
        float cpu_rem=vm1->cpu-cpu_req;
        if(cpu_rem<best_cpu){
            best_cpu=cpu_rem;
            best_vm=vm1;
        }
    }
    cnptr newCn=create_cn(cn_ogi);
    addCnInVm(newCn,best_vm);
    //printf("\nIN HEURISTICS 9");
    absent_count--;
    heuristics(t1,cn_absent,absent_count);
}

void crossover(popptr p1){
    geneptr parent1=binary_trnment(p1);
    geneptr parent2=binary_trnment(p1);

    while(parent1==parent2){
        parent2=binary_trnment(p1);
    }

    printf("\n\nparent 1 utilization:%f,parent 2 utilization:%f",parent1->utlz,parent2->utlz);

    geneptr child=(geneptr)malloc(sizeof(struct gene));
    child->hpm=NULL;
    child->next=NULL;
    child->npm=0;
    child->utlz=0;

    int min_pm=parent1->npm>parent2->npm?parent2->npm:parent1->npm;
    fitness(parent1);
    fitness(parent2);
    while(min_pm){

        //printf("\nin crsover 1");

        pmptr pm1=search_best_pm(parent1);
        pmptr pm2=search_best_pm(parent2);

        pmptr best_pm;
        geneptr best_parent,othr_parent;


        //select best pm with a probability
        int dc=(rand()%10)<PROB_PM_TO_ADD?1:0;
        if(pm1->utlz>pm2->utlz){
            if(dc){
                best_pm=pm1;

            }
            else
                best_pm=pm2;

        }
        else{
            if(dc){
                best_pm=pm2;
            }
            else{
                best_pm=pm1;
            }
        }

        if(best_pm==pm1){
            best_parent=parent1;
        }
        else{
            best_parent=parent2;
        }

        if(best_parent==parent1){
            othr_parent=parent2;
        }
        else{
            othr_parent=parent1;
        }

        //printf("\nin crsover 3");
        //copy pm into the child.
        pmptr cpm=create_pm(child);
        copy_pm(best_pm,cpm);
        //printf("\nin crsover 4");
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
        //printf("\nin crsover 5");
        //remove pm at respective parent.

        remove_pm(best_pm,best_parent);
        //printf("\nin crsover 6");
        //calculate fitness of both parents.
        fitness(parent1);
        fitness(parent2);
        min_pm--;

    }
    //printstats(child);
    fitness(child);


    //create a present containers list.
    int cn_present[CNT_SIZE];
    count_containers(child);
    cnInGene(child,cn_present);

    //unpack
    //printf("\nin crsover 6.5");
    unpack_pm(child,cn_present);
    //printf("\nin crsover 6.6");

    //create an absent containers list.
    int cn_absent[CNT_SIZE];
    int absent_count=0;
    for(int i=0;i<CNT_SIZE;i++){
        if(cn_present[i]==0){
            cn_absent[absent_count]=i;
            absent_count++;
        }
    }
    //printf("\nabsent count:%d",absent_count);

    //check the integrity of the list.1.no duplicates checked
    //remove n random containers from child.//do this in mutation.


    //add the absent containers acc to heuristic.
    //printstats(child);
    if(absent_count!=0){
        heuristics(child,cn_absent,absent_count);
        count_containers(child);
        integrity_duplicate_containers(child);
        integrity_all_containers(child);
    }
    //merge vm if less containers present.
    merge_vm(child);
    //remove empty vms
    //integrity check.1.child contains all containers.2.VM/PM size matches with the respective containers.3.no.of containers/vm/pm actually matches n

    //remove parents from genepool
    remove_gene(p1,parent1);
    //printf("\nin crsover 8");
    remove_gene(p1,parent2);
    //printf("\nin crsover 9");
    //add child to genepool.

    printstats(child);
    fitness(child);
    add_gene(p1,child);


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
void fitness(geneptr t1){
    t1->utlz=0;
    pmptr pm1=t1->hpm;
    int i=0;
    printf("\n\nCalculating fitness...");
    while(pm1!=NULL){
        vmptr vm1=pm1->hvm;
        float utlz=0;
        while(vm1!=NULL){
            cnptr cn1=vm1->hcn;
            while(cn1!=NULL){
                utlz+=cn1->cpu;
                cn1=cn1->next;
            }
            int ogi=vm1->ogi;
            float vm_overhead=data1.VM_TYPE[ogi][0]*0.1;
            utlz+=vm_overhead;
            vm1=vm1->next;
        }
        utlz=utlz/data1.pm[0];
        pm1->utlz=utlz;
        t1->utlz+=utlz;
        printf("\n\tPM:%d Utlz:%f",i,pm1->utlz);
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
            printf("\n\tVM:%d : CPU:%f/%f Mem:%f/%f No of Containers:%d",j,vm->cpu,data1.VM_TYPE[vm->ogi][0],vm->mem,data1.VM_TYPE[vm->ogi][1],vm->ncn);
            cnptr cn=vm->hcn;
            int k=0;
            while(cn!=NULL){
                printf("\n\t\tContainer:%d : CPU taken:%f Mem taken:%f",k,cn->cpu,cn->mem);
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

    float best_mem=99999;
    int vmi=-1;
    for(int i=0;i<VM_SIZE;i++){
        int vm_ogi=i;
        float vm_cpu=data1.VM_TYPE[vm_ogi][0];
        float vm_mem=data1.VM_TYPE[vm_ogi][1];

        if((data1.pm[0]<(vm_cpu*1.1)) || (data1.pm[1]<(vm_mem+200)) || (vm_mem<mem_req) || (vm_cpu<cpu_req)){
            continue;
        }
        float mem_diff=vm_mem-mem_req;
        if(best_mem>mem_diff){
            best_mem=mem_diff;
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
