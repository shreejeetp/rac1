#include<stdio.h>
#include<stdlib.h>

void rac(float containers_synthetic[200][2],float containers_real[200][2],float vm_synthetic[10][2],float vm_real[20][2],float pm[3]);

void first_fit(int pm_vm[100][100],int vm_containers[100][100],float pm[3],float vm[10][2],float containers[200][2],int testNo[3],int *testContainers);
int findFirstFitVm(float current_container[2],float pmstats[100][3],float vmstats[100][4],int pmcounter,int pm_vm[100][100],int *selectedVm);
void allocateVm(int counter[3],float current_container[3],float pmstats[100][3],float vmstats[100][4],int vm_containers[100][100],int *selectedVm);
int selectVm(float current_container[3],float vm[10][2],int *selectedVm);
void createVm(int counter[3],int pm_vm[100][100],float vm[10][2],float current_container[3],float pmstats[100][3],float vmstats[100][4],int selectedVm);
void createPm(float pmstats[100][3],float pm[3],int counter[3]);

void readcsv(
    float cont_syn[200][2],
    float cont_rea[200][2],
    float vm_syn[10][2],
    float vm_rea[20][2],
    float pm[3]
   );


int main (){
    float cont_syn[200][2];
    float cont_rea[200][2];
    float vm_syn[10][2];
    float vm_rea[20][2];
    float pm[3];
    readcsv(cont_syn,cont_rea,vm_syn,vm_rea,pm);
    
    rac(cont_syn,cont_rea,vm_syn,vm_rea,pm);

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

void rac(float containers_synthetic[200][2],float containers_real[200][2],float vm_synthetic[10][2],float vm_real[20][2],float pm[3])
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

    /*int testCn=100;
    int testVn=100;
    int testpn=100;

    int pmvm[testpn][testVn];
    int vmco[testVn][testCn][2];
    int nc[200];
    int i=0;
    int nc1[testCn];
    while(i<=testCn){
        int ran=(rand()%200);
        
        while(nc[ran]!=0){
            ran=(rand()%200);
        }
        nc[ran]=1;
        nc1[i]=ran;
        i++;
        //printf("%d, %d\n",ran,t);
    }
    
    //alloc
    int pms[3];
    pms[0]=-1;
    pms[1]=0;
    pms[2]=0;
    int vmns[3]
    vmns[0]=-1;
    vmns[1]=0;
    vmns[2]=0;
    int exn[3];
    exn[0]=0;
    int vmt=0;
    int pmvmn=0;
    for(int i=0;i<testCn;i++){
        int con_cpu=cons[nc1[i]][0];
        int con_mem=cons[nc1[i]][1];
        
        if((vmns[2]>con_mem)&&(vmns[1]>con_cpu)&&(pms[2]>con_mem)&&(pms[1]>con_cpu)){
            vmns[2]-=con_mem;
            vmns[1]-=con_cpu;
            vmco[vmns[0]][exn[0]][0]=nc1[i];
            vmco[vmns[0]][exn[0]][1]=vmt;
            exn[0]++;
        }
        else if((pms[2]>con_mem)&&(pm[1]>con_cpu)){
            //new vm
            vmns[0]++;
            exn[0]=0;
            vmt=rand()%10;
            vmns[2]=vms[vmt][1];
            vmns[1]=vms[vmt][0];
            while((vmns[2]<con_mem)&&(vmns[1]<con_cpu)){
                vmt=rand()%10;
                vmns[2]=vms[vmt][1];
                vmns[1]=vms[vmt][0];
            }
            vmns[2]-=con_mem;
            vmns[1]-=con_cpu;
            pms[2]-=con_mem;
            pm[1]-=con_cpu;
            vmco[vmns[0]][exn[0]][0]=nc1[i];
            vmco[vmns[0]][exn[0]][0]=vmt;
            pmvm[pm[0]][pmvmn++]=vmns[0];
            exn[0]++;
        }
        else{
            //new pm
            pms[0]++;
            pms[1]=pm[0];
            pms[2]=pm[1];
            //new vm
            vmns[0]++;
            exn[0]=0;
            vmt=rand()%10;
            vmns[2]=vms[vmt][1];
            vmns[1]=vms[vmt][0];
            while((vmns[2]<con_mem)&&(vmns[1]<con_cpu)){
                vmt=rand()%10;
                vmns[2]=vms[vmt][1];
                vmns[1]=vms[vmt][0];
            }
            vmns[2]-=con_mem;
            vmns[1]-=con_cpu;
            pms[2]-=con_mem;
            pms-=con_cpu;
            vmco[vmns[0]][exn[0]][0]=nc1[i];
            vmco[vmns[0]][exn[0]][0]=vmt;
            pmvm[pm[0]][pmvmn++]=vmns[0];
            exn[0]++;
        }
    }
    printf("\nNo.of PM:%d\nNo.of VM:%d\n",pmvmn,vmns[0]);
    */
