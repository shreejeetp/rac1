#include<stdio.h>
#include<stdlib.h>

//void rac(float containers_synthetic[200][2],float containers_real[200][2],float vm_synthetic[10][2],float vm_real[20][2],float pm[3]);

//void first_fit(int pm_vm[100][100],int vm_containers[100][100],float pm[3],float vm[10][2],float containers[200][2],int testNo[3],int *testContainers);
//int findFirstFitVm(float current_container[2],float pmstats[100][3],float vmstats[100][4],int pmcounter,int pm_vm[100][100],int *selectedVm);
//void allocateVm(int counter[3],float current_container[3],float pmstats[100][3],float vmstats[100][4],int vm_containers[100][100],int *selectedVm);
//int selectVm(float current_container[3],float vm[10][2],int *selectedVm);
//void createVm(int counter[3],int pm_vm[100][100],float vm[10][2],float current_container[3],float pmstats[100][3],float vmstats[100][4],int selectedVm);
//void createPm(float pmstats[100][3],float pm[3],int counter[3]);

void readcsv(
    float cont_syn[200][2],
    float cont_rea[200][2],
    float vm_syn[10][2],
    float vm_rea[20][2],
    float pm[3]
   );

void generate_random_containers(
    int n,
    int test_containers[]
   );

void first_fit
(
float pmstats[100][3],
float vmstats[100][4],
float cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float cont[200][2],
float vm[20][2],
float pm[3],
int NO_OF_TEST_CONTAINERS,
int test_containers[NO_OF_TEST_CONTAINERS]);

void printstats(
float pmstats[100][3],
float vmstats[100][4],
float cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float container[200][2]
);
void best_fit(
float pmstats[100][3],
float vmstats[100][4],
float cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float cont[200][2],
float vm[20][2],
float pm[3],
int NO_OF_TEST_CONTAINERS,
int test_containers[]
);


struct data{
    float cont_syn[200][2];
    float cont_rea[200][2];
    float vm_syn[10][2];
    float vm_rea[20][2];
    float pm[3];
};
struct population{
    float vmstats[100][4];
    float pmstats[100][3];
    float cnstats[100];
    int pmvm[100][2];
    int vmcn[100][2];
    int id_counter[5];
    int test_containers[100];
    struct population *next;
}*head_population;
typedef struct population* popptr;

void create_population(struct data data1);
void fitness(struct data data1, popptr pop);


int main (){
    struct data data1;
    //stats pm=spaceavailable,no.of vms
    //pmvm=id,vid,
    //#vm=spaceavailable,oid,no of containersa
    //#vmcn=vid,cid,
    //#cn=cid,oid
    //#idcounter=pm,vm,cn,vmcn,pmvm


    /*float vmstats[100][4];
    float pmstats[100][3];
    float cnstats[100];
    int pmvm[100][2];
    int vmcn[100][2];
    int id_counter[]={0,0,0,0,0};
    */

    readcsv(data1.cont_syn,data1.cont_rea,data1.vm_syn,data1.vm_rea,data1.pm);

    create_population(data1);

    }
void create_population(struct data data1){
    //int NO_OF_TEST_CONTAINERS=100;

    head_population=(popptr)malloc(sizeof(struct population));

    popptr t1=head_population;
    for(int i=0;i<2;i++){
        popptr t2=(popptr)malloc(sizeof(struct population));
        t1->next=t2;
        generate_random_containers(100,t2->test_containers);

        first_fit(t2->pmstats,t2->vmstats,t2->cnstats,t2->pmvm,t2->vmcn,t2->id_counter,data1.cont_syn,data1.vm_syn,data1.pm,100,t2->test_containers);

        fitness(data1,t2);
        printstats(t2->pmstats,t2->vmstats,t2->cnstats,t2->pmvm,t2->vmcn,t2->id_counter,data1.cont_syn);
        t2->next=NULL;
        t1=t2;
    }
}

void fitness(struct data data1, popptr pop){
    float energy[pop->id_counter[0]];
    for(int i=0;i<pop->id_counter[0];i++){
            energy[i]=0;
            int pid=i;
            for(int j=0;j<pop->id_counter[4];j++){
                if(pop->pmvm[j][0]==pid){
                    int vid=pop->pmvm[j][1];
                    int vm_oid=pop->vmstats[j][2];
                    float vm_overhead=data1.vm_syn[vm_oid][0]*.1;
                    float cont_cpu=0;
                    for(int k=0;k<pop->id_counter[3];k++){
                        if(pop->vmcn[k][0]==vid){
                            int cid=pop->vmcn[k][1];
                            int c_oid=pop->cnstats[cid];
                            cont_cpu+=data1.cont_syn[c_oid][0];
                        }
                    }
                    energy[i]+=vm_overhead+cont_cpu;
                }
            }
            energy[i]=energy[i]/data1.pm[0];
            printf("\nutilization of pm %d is %f",i,energy[i]);
    }
}

void generate_random_containers(
    int n,
    int test_containers[]
   ){
    int i=0;
    static int check_bucket[200];
    while(i<n){
        int ran=(rand()%200);
        while(check_bucket[ran]!=0){
            ran=(rand()%200);
        }
        test_containers[i]=ran;
        check_bucket[ran]=1;
        i++;
    }
   }
void readcsv(
    float cont_syn[200][2],
    float cont_rea[200][2],
    float vm_syn[10][2],
    float vm_rea[20][2],
    float pm[3]
   ){

    FILE *fpr;
    FILE *fpr2;
    FILE *fpr3;
    FILE *fpr4;
    FILE *fpr5;
    char c;
    fpr=fopen("data/Container200_ten.csv","r");
    fpr2=fopen("data/Container200_twenty.csv","r");
    fpr3=fopen("data/VMConfig_ten.csv","r");
    fpr4=fopen("data/VMConfig_twenty.csv","r");
    fpr5=fopen("data/PMConfig_small.csv","r");

    if(fpr==NULL||fpr2==NULL||fpr3==NULL||fpr4==NULL||fpr5==NULL){
        printf("\n ERROR in opening file\n");
    }
    //reading container file for synthetic vms
    int i=0;
    while(i<200){
        fscanf(fpr,"%f,%f\n",&cont_syn[i][0],&cont_syn[i][1]);
        i++;
    }
    //reading container file for real vms
    i=0;
    while(i<200){
        fscanf(fpr2,"%f,%f\n",&cont_rea[i][0],&cont_rea[i][1]);
        i++;
    }

    //reading synthetic vm config
    i=0;
    while(i<10){
        fscanf(fpr3,"%f,%f\n",&vm_syn[i][0],&vm_syn[i][1]);
        i++;
    }

    //reading real vm config
    i=0;
    while(i<20){
        fscanf(fpr4,"%f,%f\n",&vm_rea[i][0],&vm_rea[i][1]);
        i++;
    }

    //reading pm config
    fscanf(fpr5,"%f\n%f\n%f",&pm[0],&pm[1],&pm[2]);
    fclose(fpr);
    fclose(fpr2);
    fclose(fpr3);
    fclose(fpr4);
    fclose(fpr5);

   }

void first_fit(
float pmstats[100][3],
float vmstats[100][4],
float cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float cont[200][2],
float vm[20][2],
float pm[3],
int NO_OF_TEST_CONTAINERS,
int test_containers[]
){

    for(int i=0;i<NO_OF_TEST_CONTAINERS;i++){
        int current_container=test_containers[i];

        //if theres space in any vm
        int selected_vid=-1;
        for(int j=0;j<id_counter[1];j++){
            float temp_cpu=vmstats[j][0]-cont[current_container][0];
            float temp_mem=vmstats[j][1]-cont[current_container][1];

            if((temp_cpu>=0)&&(temp_mem>=0)){
                selected_vid=j;
                break;
            }
        }
        if(selected_vid!=-1){
            cnstats[id_counter[2]]=current_container;
            int cid=id_counter[2];
            int vmcnid=id_counter[3];

            vmstats[selected_vid][0]-=cont[current_container][0];
            vmstats[selected_vid][1]-=cont[current_container][1];
            vmstats[selected_vid][3]+=1;

            vmcn[vmcnid][0]=selected_vid;
            vmcn[vmcnid][1]=cid;


            id_counter[2]++;
            id_counter[3]++;
            }
        else{
            //select random vm
            int selected_vm=(rand()%20);
            while(( cont[current_container][0] > vm[selected_vm][0] ) || ( cont[current_container][1] > vm[selected_vm][1] ) || ( vm[selected_vm][0] * 1.1 > pm[0] ) || (( vm[selected_vm][1] + 200 ) > pm[1] )){
                selected_vm=(rand()%20);
            }
            //see if there's space in any existing pm
            int selected_pid=-1;
            for(int pid=0;pid<id_counter[0];pid++){
                if(( pmstats[pid][0] >= vm[selected_vm][0]*1.1 ) && ( pmstats[pid][1] >= vm[selected_vm][1] + 200 )){
                    selected_pid=pid;
                    break;
                }
            }
            int pid=selected_pid;

            if(pid==-1){
                pid=id_counter[0];
                id_counter[0]+=1;

                pmstats[pid][0]=pm[0];
                pmstats[pid][1]=pm[1];
                pmstats[pid][2]=0;

            }
            int vid=id_counter[1];
            id_counter[1]++;
            vmstats[vid][0]=vm[selected_vm][0];
            vmstats[vid][1]=vm[selected_vm][1];
            vmstats[vid][2]=selected_vm;
            vmstats[vid][3]=0;

            pmstats[pid][0]=pmstats[pid][0]-(vm[selected_vm][0]*1.1);
            pmstats[pid][1]-=vm[selected_vm][1]+200;
            pmstats[pid][2]+=1;
            pmvm[id_counter[4]][0]=pid;
            pmvm[id_counter[4]][1]=vid;

            id_counter[4]++;

            int cid=id_counter[2];

            vmstats[vid][0]-=cont[current_container][0];
            vmstats[vid][1]-=cont[current_container][1];
            vmstats[vid][3]+=1;

            vmcn[id_counter[3]][0]=vid;
            vmcn[id_counter[3]][1]=cid;


            id_counter[2]++;
            id_counter[3]++;
        }
    }
}

void best_fit(
float pmstats[100][3],
float vmstats[100][4],
float cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float cont[200][2],
float vm[20][2],
float pm[3],
int NO_OF_TEST_CONTAINERS,
int test_containers[]
){

    for(int i=0;i<NO_OF_TEST_CONTAINERS;i++){
        int current_container=test_containers[i];

        //if theres space in any vm
        int selected_vid=-1;
        float mincpu=99999;
        for(int j=0;j<id_counter[1];j++){
            float temp_cpu=vmstats[j][0]-cont[current_container][0];
            float temp_mem=vmstats[j][1]-cont[current_container][1];

            if((temp_cpu>=0)&&(temp_mem>=0)&&(mincpu>temp_cpu)){
                selected_vid=j;
                mincpu=temp_cpu;
            }
        }
        if(selected_vid!=-1){
            cnstats[id_counter[2]]=current_container;
            int cid=id_counter[2];
            int vmcnid=id_counter[3];

            vmstats[selected_vid][0]-=cont[current_container][0];
            vmstats[selected_vid][1]-=cont[current_container][1];
            vmstats[selected_vid][3]+=1;

            vmcn[vmcnid][0]=selected_vid;
            vmcn[vmcnid][1]=cid;


            id_counter[2]++;
            id_counter[3]++;
            }
        else{
            //select random vm
            int selected_vm=(rand()%20);
            while(( cont[current_container][0] > vm[selected_vm][0] ) || ( cont[current_container][1] > vm[selected_vm][1] ) || ( vm[selected_vm][0] * 1.1 > pm[0] ) || (( vm[selected_vm][1] + 200 ) > pm[1] )){
                selected_vm=(rand()%20);
            }
            //see if there's space in any existing pm
            int selected_pid=-1;
            for(int pid=0;pid<id_counter[0];pid++){
                if(( pmstats[pid][0] >= vm[selected_vm][0]*1.1 ) && ( pmstats[pid][1] >= vm[selected_vm][1] + 200 )){
                    selected_pid=pid;
                    break;
                }
            }
            int pid=selected_pid;

            if(pid==-1){
                pid=id_counter[0];
                id_counter[0]+=1;

                pmstats[pid][0]=pm[0];
                pmstats[pid][1]=pm[1];
                pmstats[pid][2]=0;

            }
            int vid=id_counter[1];
            id_counter[1]++;
            vmstats[vid][0]=vm[selected_vm][0];
            vmstats[vid][1]=vm[selected_vm][1];
            vmstats[vid][2]=selected_vm;
            vmstats[vid][3]=0;

            pmstats[pid][0]=pmstats[pid][0]-(vm[selected_vm][0]*1.1);
            pmstats[pid][1]-=vm[selected_vm][1]+200;
            pmstats[pid][2]+=1;
            pmvm[id_counter[4]][0]=pid;
            pmvm[id_counter[4]][1]=vid;

            id_counter[4]++;

            int cid=id_counter[2];

            vmstats[vid][0]-=cont[current_container][0];
            vmstats[vid][1]-=cont[current_container][1];
            vmstats[vid][3]+=1;

            vmcn[id_counter[3]][0]=vid;
            vmcn[id_counter[3]][1]=cid;


            id_counter[2]++;
            id_counter[3]++;
        }
    }
}

void printstats(
float pmstats[100][3],
float vmstats[100][4],
float cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float container[200][2]
){
    //for(int i=0;i<id_counter[3];)
    for(int i=0;i<id_counter[0];i++){
        printf("\n PM:%d Cpu left:%f Mem left:%f No. of VMs:%d VM:\n",i,pmstats[i][0],pmstats[i][1],pmstats[i][2]);
        for(int j=0;j<id_counter[4];j++){
            if(pmvm[j][0]==i){
                printf("\n\t VM:%d Cpu left:%f Mem left:%f No. of Containers:%d Containers:\n",pmvm[j][1],vmstats[pmvm[j][1]][0],vmstats[pmvm[j][1]][1],vmstats[pmvm[j][1]][3]);
                for(int k=0;k<id_counter[3];k++){
                    if(vmcn[k][0]==pmvm[j][1]){
                        int cont=(int)cnstats[vmcn[k][1]];
                        printf("\n\t\t Container:%d Cpu used:%f Mem used:%f.",vmcn[k][1],container[cont][0],container[cont][1]);
                    }
                }
            }
        }
    }
}
/*void rac(float containers_synthetic[200][2],float containers_real[200][2],float vm_synthetic[10][2],float vm_real[20][2],float pm[3])
{
    //declaring test variables
    int testNo[3]={100,100,100};
    //0-container,1-vm,2-pm
    int testContainers[testNo[0]];

     //generating random numbers
    int i=0;
    while(i<testNo[0]){
        int ran=(rand()%200);
        testContainers[i]=ran;
        i++;
    }
    int pm_vm[testNo[2]][testNo[1]];
    int vm_containers[testNo[1]][testNo[0]];

    first_fit(pm_vm,vm_containers,pm,vm_synthetic,containers_synthetic,testNo,testContainers);


}




void first_fit(int pm_vm[100][100],int vm_containers[100][100],float pm[3],float vm[10][2],float containers[200][2],int testNo[3],int *testContainers){

    //0-cpu
    //1-mem
    //2-index to vm[]//type of vm
    //3-no. of containers
    float vmstats[100][4];

    //0,1-same
    //2-no.of vms inside
    float pmstats[100][3];

    //0-container //should be reset on new vm
    //1-vm      //vm%pm=vm counter for pm
    //2-pm
    //allocating first container
    int counter[3];
    counter[0]=0;
    counter[1]=0;
    counter[2]=0;


    pmstats[0][0]=pm[0];
    pmstats[0][1]=pm[1];
    pmstats[0][2]=1;

    int ranvm=rand()%10;
    vmstats[0][0]=vm[ranvm][0];
    vmstats[0][1]=vm[ranvm][1];
    vmstats[0][2]=ranvm;
    vmstats[0][3]=1.00;

    pm_vm[counter[2]][counter[1]]=0;

    vm_containers[0][0]=*(testContainers);

    for(int i=1;i<testNo[0];i++){
        float current_container[3];
        current_container[0]=containers[*(testContainers+i)][0];
        current_container[1]=containers[*(testContainers+i)][1];
        current_container[2]=*(testContainers+i);

        //printf("\n%f\n",pmstats[counter[2]][0]);

        //if current pm has space
        if((pmstats[counter[2]][0]>current_container[0])&&(pmstats[counter[2]][1]>current_container[1])){
            //find a fit vm
            int selectedVm;
            int found=findFirstFitVm(current_container,pmstats,vmstats,counter[2],pm_vm,&selectedVm);
            //allocate if found
            if(found){
                    allocateVm(counter,current_container,pmstats,vmstats,vm_containers,&selectedVm);
                    //printf("\n%d %d %d\n",counter[0],counter[1],counter[2]);
            }
            else{
                //select a suitable random new vm
                int selectedVm;
                int found1=selectVm(current_container,vm,&selectedVm);
                if(found1){
                //if current pm can bear overhead
                    if((pmstats[counter[2]][0]>(vm[selectedVm][0]*1.1))&&(pmstats[counter[2]][1]>(vm[selectedVm][1]+200))){
                        //allocatePm(pmstats,vm,pm,selectedVm);
                        createVm(counter,pm_vm,vm,current_container,pmstats,vmstats,selectedVm);
                        selectedVm=counter[1];
                        allocateVm(counter,current_container,pmstats,vmstats,vm_containers,&selectedVm);
                    }
                    //create new pm
                    else{
                        int selectedVm;
                        createPm(pmstats,pm,counter);
                        createVm(counter,pm_vm,vm,current_container,pmstats,vmstats,selectedVm);
                        selectedVm=counter[1];
                        allocateVm(counter,current_container,pmstats,vmstats,vm_containers,&selectedVm);
                    }
                }
                else{
                    printf("\nNo Vm with space required available\n");
                    exit(0);
                }
            }
        }
        //create new pm
        else{
            int selectedVm;
            createPm(pmstats,pm,counter);
            selectVm(current_container,vm,&selectedVm);
            createVm(counter,pm_vm,vm,current_container,pmstats,vmstats,selectedVm);
            allocateVm(counter,current_container,pmstats,vmstats,vm_containers,&selectedVm);
        }
    }
    //printf("\n%f\n",pmstats[counter[2]][0]);

    for(int i=0;i<=counter[2];i++){
        printf("\n%d PM(%f cpu remaining, %f mem remaining) of (%f and %f):\n",i+1,pmstats[i][0],pmstats[i][1],pm[0],pm[1]);
        for(int j=0;j<pmstats[i][2];j++){
            int vmn=(int)vmstats[pm_vm[i][j]][2];
            printf("\n%d VM (%f cpu remaining, %f mem remaining) of (%f and %f) :\n",j+1,vmstats[pm_vm[i][j]][0],vmstats[pm_vm[i][j]][1],vm[vmn][0],vm[vmn][1]);
            for(int k=0;k<=vmstats[j][3];k++){
                printf("%d Container : %f cpu %f mem\n",k+1,containers[vm_containers[j][k]][0],containers[vm_containers[j][k]][1]);
            }
        }

    }
}

int findFirstFitVm(float current_container[2],float pmstats[100][3],float vmstats[100][4],int pmcounter,int pm_vm[100][100],int *selectedVm){
    int n=(int)pmstats[pmcounter][2];
    int i=0;
    while(i<n){
        if((vmstats[pm_vm[pmcounter][i]][0]>=current_container[0])&&(vmstats[pm_vm[pmcounter][i]][1]>=current_container[1])){
            *selectedVm=pm_vm[pmcounter][i];
            //*selectedVm=i;
            return(1);
        }
    i++;
    }
    return(0);
}

void allocateVm(int counter[3],float current_container[3],float pmstats[100][3],float vmstats[100][4],int vm_containers[100][100],int *selectedVm){
    counter[0]++;
    vm_containers[*selectedVm][(int)vmstats[*selectedVm][3]]=current_container[2];
    //decrementing the cpu and mem available
    vmstats[*selectedVm][0]-=current_container[0];
    vmstats[*selectedVm][1]-=current_container[1];
    vmstats[*selectedVm][3]++;
    pmstats[counter[2]][0]-=(current_container[0]);
    pmstats[counter[2]][1]-=(current_container[1]);

}

void createVm(int counter[3],int pm_vm[100][100],float vm[10][2],float current_container[3],float pmstats[100][3],float vmstats[100][4],int selectedVm){
    counter[1]++;
    counter[0]=-1;
    pmstats[counter[2]][2]++;
    pm_vm[counter[2]][counter[1]%(int)pmstats[counter[2]][2]]=counter[1];
    vmstats[counter[1]][0]=vm[selectedVm][0];
    vmstats[counter[1]][1]=vm[selectedVm][1];
    vmstats[counter[1]][2]=selectedVm;
    vmstats[counter[1]][3]=0;
    pmstats[counter[2]][0]-=vm[selectedVm][0]*0.1;
    pmstats[counter[2]][1]-=200;
}

int selectVm(float current_container[3],float vm[10][2],int *selectedVm){
    int i=rand()%10;
    while((vm[i][0]<current_container[0])||(vm[i][1]<current_container[1])){
        i=rand()%10;
    }
    *selectedVm=i;
    return 1;
}

void createPm(float pmstats[100][3],float pm[3],int counter[3]){
    counter[2]++;
    pmstats[counter[2]][0]=pm[0];
    pmstats[counter[2]][1]=pm[1];
    pmstats[counter[2]][2]=0;

}

 */
