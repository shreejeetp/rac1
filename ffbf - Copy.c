
/*
Todo
    -Each time generate Randomize test_containers positions out of selected containers. so data1 will have test containers and each popptr will have their own random positon version of it.
instead of PM utilization and allocation, allocate VM.
    +Read all the code, ensure that there is no corelation between id_counter[2] and number of containers in the machine.
    +remove 3rd index in id_counter, redundant
    +make dynamic sizes of population data.
    +increase containers samples so that more pms would be created.
    +create a function that creates a containers array from combination of all different real and synthetic containers.
    +merge lower utilized pms and vms in child
    +remove empty vms.
    +make seperate functions to remove a vm, pm, and a container to maintain integrity.
Bugs
    -In crossover, it only copies containers of one gene pool.
    +In create population, last gene isnt counter as next of second last is set to NULL somehow
    After remove_containers(),the containers are removed from between and the numbering might cause problem somewhere.
*/


#include<stdio.h>
#include<stdlib.h>
#define NO_POP 10
#define NO_CROSSOVER 3
//void rac(float containers_synthetic[200][2],float containers_real[200][2],float vm_synthetic[10][2],float vm_real[20][2],float pm[3]);

//void first_fit(int pm_vm[100][100],int vm_containers[100][100],float pm[3],float vm[10][2],float containers[200][2],int testNo[3],int *testContainers);
//int findFirstFitVm(float current_container[2],float pmstats[100][3],float vmstats[100][4],int pmcounter,int pm_vm[100][100],int *selectedVm);
//void allocateVm(int counter[3],float current_container[3],float pmstats[100][3],float vmstats[100][4],int vm_containers[100][100],int *selectedVm);
//int selectVm(float current_container[3],float vm[10][2],int *selectedVm);
//void createVm(int counter[3],int pm_vm[100][100],float vm[10][2],float current_container[3],float pmstats[100][3],float vmstats[100][4],int selectedVm);
//void createPm(float pmstats[100][3],float pm[3],int counter[3]);

void readcsv();

void generate_random_containers();
void generate_random_containers2(
    int test_containers[]
);

void first_fit
(
float pmstats[100][4],
float vmstats[100][4],
int cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float cont[200][2],
float vm[10][2],
float pm[3],
int NO_OF_TEST_CONTAINERS,
int test_containers[NO_OF_TEST_CONTAINERS]);

void printstats(
float pmstats[100][4],
float vmstats[100][4],
int cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float container[200][2]
);

void best_fit(
float pmstats[100][4],
float vmstats[100][4],
int cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float cont[200][2],
float vm[20][2],
float pm[3],
int NO_OF_TEST_CONTAINERS,
int test_containers[]
);
void search_best_gene();

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
struct population{
    float vmstats[200][4];//0-cpu remaining,1-mem remaining,2-og vm index,3-no. of containers in it
    float pmstats[200][4];//0-cpu remaining,1-mem remaining,2-no.of vm in it,3-utlz
    int cnstats[200];//og test_container index
    int pmvm[200][2];//0-pm,1-vm
    int vmcn[200][2];//0-vm,1-cn
    int id_counter[5];//0-pm,1-vm,2-cont,3-vmcn,4-pmvm
    float utlz;
    int test_containers[100];
    int total;
    struct population *next;
}*head_population;
typedef struct population* popptr;


void merge_vm(popptr pop);
void remove_vm(popptr pop, int vm);
void create_population(struct data data1);
void fitness(popptr pop);
void remove_containers(popptr pop,int bucket[100]);
void remove_container(popptr pop,int sc);
void replace_container(popptr pop,int sc,int sc1,int vm1);
void remove_gene(popptr pop);
void add_gene(popptr pop);
popptr binary_selection();
void printstats1(popptr t2);
void crossover();
void heuristics(popptr child,int left_out[100],int left_out_count);


int main (){
    //stats pm=spaceavailable,no.of vms
    //pmvm=id,vid,
    //#vm=spaceavailable,oid,no of containers
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

    readcsv();

    create_population(data1);
    printf("\ntotal %d genes",head_population->total);
    for(int i=0;i<NO_CROSSOVER;i++){
        crossover();
    }
    search_best_gene();
    printf("\ntotal %d genes",head_population->total);
    printf("\nend");
}
void crossover(){
    popptr parent1=binary_selection();
    //printstats1(selected);
    popptr parent2=binary_selection();
    while(parent1->utlz==parent2->utlz){
        parent2=binary_selection();
    }
    printf("\n\nparent 1 utilization:%f,parent 2 utilization:%f",parent1->utlz,parent2->utlz);

    popptr child=(popptr)malloc(sizeof(struct population));

    /*child->id_counter[0]=0;
    child->id_counter[1]=0;
    child->id_counter[2]=0;
    child->id_counter[3]=0;
    child->id_counter[4]=0;
    */

    fitness(parent1);
    fitness(parent2);
    int min_pm=parent1->id_counter[0]<parent2->id_counter[0]?parent1->id_counter[0]:parent2->id_counter[0];
    int removal_bucket[100];
    for(int i=0;i<100;i++){
        removal_bucket[i]=0;
    }
    while(min_pm){
        //find best utilized pm in both parents
        //0,x:parent selector,x,0:pm selector,x,1:utilization

        float best_utlz=0;
        int best_utlz_pm=0;
        for(int i=0;i<parent1->id_counter[0];i++){
            if(parent1->pmstats[i][3]>best_utlz){
                best_utlz=parent1->pmstats[i][3];
                best_utlz_pm=i;
            }
        }
        int good_parent=0;
        for(int i=0;i<parent2->id_counter[0];i++){
            if(parent2->pmstats[i][3]>best_utlz){
                best_utlz=parent2->pmstats[i][3];
                best_utlz_pm=i;
                good_parent=1;
            }
        }


        if(good_parent==1){
            popptr temp=parent1;
            parent1=parent2;
            parent2=temp;
        }
        printf("\nPm with utlz :%f is selected\n",best_utlz);
        //add it to the child
        int container_bucket[100];
        for(int i=0;i<100;i++){
            container_bucket[i]=0;
        }

        child->pmstats[child->id_counter[0]][0]=parent1->pmstats[best_utlz_pm][0];
        child->pmstats[child->id_counter[0]][1]=parent1->pmstats[best_utlz_pm][1];
        child->pmstats[child->id_counter[0]][2]=parent1->pmstats[best_utlz_pm][2];
        child->pmstats[child->id_counter[0]][3]=parent1->pmstats[best_utlz_pm][3];

        for(int i=0;i<parent1->id_counter[4];i++){
            if(parent1->pmvm[i][0]==best_utlz_pm){

                if(parent1->pmvm[i][1]==-1){continue;}
                int selected_vm=parent1->pmvm[i][1];


                child->pmvm[child->id_counter[4]][0]=child->id_counter[0];
                child->pmvm[child->id_counter[4]][1]=child->id_counter[1];
                child->id_counter[4]++;


                child->vmstats[child->id_counter[1]][0]=parent1->vmstats[selected_vm][0];
                child->vmstats[child->id_counter[1]][1]=parent1->vmstats[selected_vm][1];
                child->vmstats[child->id_counter[1]][2]=parent1->vmstats[selected_vm][2];
                child->vmstats[child->id_counter[1]][3]=parent1->vmstats[selected_vm][3];


                for(int j=0;j<parent1->id_counter[2];j++){

                    if((parent1->vmcn[j][0]==selected_vm)&&(parent1->vmcn[j][1]!=-1)){
                        int selected_cn=parent1->vmcn[j][1];
                        child->vmcn[child->id_counter[3]][0]=child->id_counter[1];
                        child->vmcn[child->id_counter[3]][1]=child->id_counter[2];

                        child->cnstats[child->id_counter[2]]=parent1->cnstats[selected_cn];
                        //printf("\n%d,%d added to the child\n",parent1->cnstats[selected_cn],selected_cn);
                        container_bucket[parent1->cnstats[selected_cn]]=1;
                        removal_bucket[parent1->cnstats[selected_cn]]=1;

                        child->id_counter[2]++;
                        child->id_counter[3]++;

                    }
                }
                child->id_counter[1]++;
            }

        }
        child->id_counter[0]++;

        //remove the pm from parent

        parent1->pmstats[best_utlz_pm][3]=0;


        remove_containers(parent1,container_bucket);

        //remove the same containers from second PM
        remove_containers(parent2,container_bucket);
        //calculate utilization stats again for both
        min_pm--;
    }

    //mutation, remove 5 random containers from child



    int removal_bucket_child[100];
    for(int i=0;i<100;i++){
        removal_bucket_child[i]=0;
    }
    for(int i=0;i<15;i++){
        int ran=rand()%100;
        while((removal_bucket_child[ran]==1)||(removal_bucket[ran]==0)){
            ran=rand()%100;
        }
        //printf("\n%d added in removal bucket\n",ran);
        removal_bucket_child[ran]=1;

    }

    fitness(child);

    remove_containers(child,removal_bucket_child);


    //printstats1(child);
    fitness(child);

    //find left out containers, diff btw test_containers and bucket here
    int left_out[100];
    int left_out_count=0;
    for(int i=0;i<100;i++){
        left_out[i]=0;
    }
    for(int i=0;i<100;i++){
        if((removal_bucket[i]==0)||(removal_bucket_child[i]==1)){
            left_out[i]=1;
            left_out_count++;
            //printf("\n%d was left out\n",i);
        }
    }
    int left_out1[100];
    int lft=0;
    for(int i=0;i<100;i++){
        if(left_out[i]){
            left_out1[lft]=i;
            lft++;
        }
    }
    printf("\ntotal %d were left out of child\n",left_out_count);

    //add left out containers back using heuristics
    heuristics(child,left_out1,left_out_count);
    fitness(child);
    //remove parent1 and 2 from populations, add child
    remove_gene(parent1);
    remove_gene(parent2);
    add_gene(child);
    merge_vm(child);
    printstats1(child);
}
void merge_vm(popptr pop){
    int MBT=5;//merge below this number
    for(int i=0;i<pop->id_counter[1];i++){
        if((pop->vmstats[i][3]>=MBT)||(pop->vmstats[i][2]==-1)){
            continue;
        }
        int cont_count=0;
        int bucket[100];
        int heu_bucket[100];
        for(int j=0;j<100;j++){
            bucket[j]=0;
            heu_bucket[j]=0;
        }
        for(int j=0;j<pop->id_counter[3];j++){
            int vm=pop->vmcn[j][0];
            int cn=pop->vmcn[j][1];
            if(vm!=i){continue;}

            bucket[pop->cnstats[cn]]=1;
            heu_bucket[cont_count]=pop->cnstats[cn];
            cont_count++;
        }
        remove_containers(pop,bucket);
        remove_vm(pop,i);
        heuristics(pop,heu_bucket,cont_count);
    }
}
void remove_vm(popptr pop, int vm){
    for(int i=0;i<pop->id_counter[4];i++){
        int vm1=pop->pmvm[i][1];
        int pm1=pop->pmvm[i][0];
        if(vm1!=vm){
            continue;
        }
        pop->pmvm[i][1]=-1;
        int ogvm=(int)pop->vmstats[vm][2];
        pop->pmstats[pm1][0]+=(data1.vm_syn[ogvm][0]*1.1);
        pop->pmstats[pm1][1]+=(data1.vm_syn[ogvm][1]+200);
        pop->pmstats[pm1][2]--;

        pop->vmstats[vm][2]=-1;

    }
}
void add_gene(popptr pop){
    popptr temp1=head_population->next;
    while(temp1->next!=NULL){
        temp1=temp1->next;
    }
    temp1->next=pop;
    pop->next=NULL;
    head_population->total++;
}
void remove_gene(popptr pop){
    popptr temp1=head_population->next;
    popptr temp2=head_population;
    while(temp1!=NULL){
        if(pop==temp1){
            temp2->next=temp1->next;
            head_population->total--;
            break;
        }
        temp2=temp1;
        temp1=temp1->next;

    }
}
void heuristics(popptr child,int left_out[100],int left_out_count){
    if(left_out_count==0){
        return;
    }
    //printf("\nallocating %d of left out",left_out_count);
    //always allocate in the pm which has higher utilization

    float best_utlz=0;
    int best_utlz_pm=0;

    for(int i=0;i<child->id_counter[0];i++){
            if(child->pmstats[i][3]>best_utlz){
                best_utlz=child->pmstats[i][3];
                best_utlz_pm=i;
            }
        }
    int top=left_out_count-1;

    float cpu_req=data1.cont_syn200[data1.test_containers[left_out[top]]][0];
    float mem_req=data1.cont_syn200[data1.test_containers[left_out[top]]][1];

    //check if this container fills any vm completely
    float complete_fill_threshhold=99.99;
    int cmplt_fill=0;
    int bttr_fit=0;

    //printf("\ntrying complete fit");
    for(int i=0;i<child->id_counter[1];i++){
        //use pmvmpmvm here later so that it is allocated only in best utlz pm
        if((cpu_req>child->vmstats[i][0])||(mem_req>child->vmstats[i][1])||(child->vmstats[i][2]==-1)){
            continue;
        }
        float cpu_rem=child->vmstats[i][0]-cpu_req;
        if(cpu_rem<=complete_fill_threshhold){
            cmplt_fill=1;
            //allocate_cont(child,best_utlz_pm,i,left_out[top]);
            int idc=child->id_counter[2];
            child->vmcn[idc][0]=i;
            child->vmcn[idc][1]=idc;

            child->vmstats[i][3]++;
            child->vmstats[i][0]-=cpu_req;
            child->vmstats[i][1]-=mem_req;
            child->cnstats[idc]=left_out[top];

            child->id_counter[2]++;
            child->id_counter[3]++;
            printf("\ncomplete fit used");
            left_out_count--;
            break;
        }
    }

    //check if this container fits better
    if(!cmplt_fill){

        for(int i=0;i<child->id_counter[1];i++){
            if(child->vmstats[i][2]==-1){continue;}
            for(int pmi=0;pmi<child->id_counter[4];pmi++){
                if((child->pmvm[pmi][1]==i)&&(child->pmvm[pmi][0]==best_utlz_pm)){

                    if(!bttr_fit){
                        for(int j=0;j<child->id_counter[2];j++){
                            //printf("\ntrying bttr fit");
                            if(child->vmcn[j][0]!=i){
                                continue;
                            }
                            int testIn=child->cnstats[child->vmcn[j][1]];
                            int contIn=data1.test_containers[testIn];

                            if((data1.cont_syn200[contIn][0]-cpu_req>=0)){
                                continue;
                            }

                            float cpu_diff=child->vmstats[i][0]+data1.cont_syn200[contIn][0]-cpu_req;
                            float mem_diff=child->vmstats[i][1]+data1.cont_syn200[contIn][1]-mem_req;

                            if((cpu_diff<0)||(mem_diff<0)){
                                continue;
                            }
                            //printf("\n%d,rep mem:%f,mem_diff:%f",i,data1.cont_syn200[contIn][1],mem_diff);
                            printf("\n%d,%d,%d",i,testIn,left_out[top]);
                            printf("\nbetter fit used");
                            bttr_fit=1;
                            replace_container(child,testIn,left_out[top],i);
                            left_out[top]=testIn;

                            break;
                        }
                    }
                }
            }
        }
    }
    //printf("\n1. pm:%d",child->id_counter[0]);
    //allocate using best fit if not
    int best_fit=!cmplt_fill&&!bttr_fit;
    //printf("\n%d,%d,%d",cmplt_fill,bttr_fit,best_fit);
    if(best_fit){
        printf("\nbest fit used");
        int best_vm=-1;
        float min_cpu=9999;
        for(int i=0;i<child->id_counter[1];i++){
            if((cpu_req>child->vmstats[i][0])||(mem_req>child->vmstats[i][1])||(child->vmstats[i][2]==-1)){
                //printf("\nbestfit:cpu %f - %f,mem %f - %f",cpu_req,child->vmstats[i][0],mem_req,child->vmstats[i][1]);
                continue;
            }
            if(child->vmstats[i][0]-cpu_req<min_cpu){
                min_cpu=child->vmstats[i][0]-cpu_req;
                best_vm=i;
            }
        }
        if(best_vm<0){
            printstats1(child);
            printf("\nerror:cant find vm in heuristics best fit,%d",data1.cont_syn200[data1.test_containers[left_out[top]]]);
            exit(0);
        }
        int idc=child->id_counter[2];

        child->vmcn[idc][0]=best_vm;
        //child->vmcn[child->id_counter[2]][0]=best_vm;

        child->vmcn[idc][1]=idc;


        child->vmstats[best_vm][3]++;
        child->vmstats[best_vm][0]-=cpu_req;
        child->vmstats[best_vm][1]-=mem_req;

        child->cnstats[idc]=left_out[top];

        child->id_counter[2]++;
        child->id_counter[3]++;

        left_out_count--;

    }
    //printf("\n4. pm:%d",child->id_counter[0]);
    heuristics(child,left_out,left_out_count);

}
void replace_container(popptr pop,int sc,int sc1,int vm1){
    int counter=0;
    for(int j=0;j<pop->id_counter[3];j++){
            int vm=pop->vmcn[j][0];
            int cn=pop->vmcn[j][1];
            if(pop->cnstats[cn]==sc){
                            //printf("\n%d removed\n",sc);
                            pop->cnstats[cn]=sc1;
                            pop->vmstats[vm][0]+=data1.cont_syn200[data1.test_containers[sc]][0];
                            pop->vmstats[vm][1]+=data1.cont_syn200[data1.test_containers[sc]][1];
                            pop->vmstats[vm][0]-=data1.cont_syn200[data1.test_containers[sc1]][0];
                            pop->vmstats[vm][1]-=data1.cont_syn200[data1.test_containers[sc1]][1];

                            if(vm1!=vm){
                                //printf("\nerror,%d, %f - %f\n",pop->vmcn[j][0],data1.cont_syn200[data1.test_containers[sc]][1], pop->vmstats[pop->vmcn[j][0]][1]);
                                printf("\nerror,%d, %d - %d\n\n\n",vm,sc,sc1);
                                //exit(0);
                            }
                            return;

        }
    }


}

void remove_containers(popptr pop,int bucket[100]){
    int counter=0;

    for(int sc=0;sc<100;sc++){
        if(bucket[sc]==1){
            for(int i=0;i<pop->id_counter[2];i++){
                if(pop->cnstats[i]==sc){
                    for(int j=0;j<pop->id_counter[3];j++){
                        if(pop->vmcn[j][1]==i){
                            //printf("\n%d removed\n",sc);
                            pop->vmcn[j][1]=-1;
                            pop->vmstats[pop->vmcn[j][0]][0]+=data1.cont_syn200[data1.test_containers[sc]][0];
                            pop->vmstats[pop->vmcn[j][0]][1]+=data1.cont_syn200[data1.test_containers[sc]][1];
                            pop->vmstats[pop->vmcn[j][0]][3]--;
                            counter++;
                        }
                    }
                }
            }
        }
    }
    printf("\n%d containers removed.\n",counter);
}
void create_population(struct data data1){
    //int NO_OF_TEST_CONTAINERS=100;

    head_population=(popptr)malloc(sizeof(struct population));
    head_population->total=0;
    popptr t1=head_population;

    generate_random_containers();

    for(int i=0;i<NO_POP;i++){
        popptr t2=(popptr)malloc(sizeof(struct population));
        t1->next=t2;
        //printf("\n\n1 %d",i);

        //printf("\n\n2 %d",i);
        generate_random_containers2(t2->test_containers);
        //printf("\n\n2.1 %d",i);
        first_fit(t2->pmstats,t2->vmstats,t2->cnstats,t2->pmvm,t2->vmcn,t2->id_counter,data1.cont_syn200,data1.vm_syn,data1.pm,100,t2->test_containers);

        //printf("\n\n3 %d",i);
        fitness(t2);
        //printf("\n\n4 %d",i);
        printstats1(t2);
        //printf("\n\n5 %d",i);
        //fflush(stdout);
        t2->next=NULL;
        t1=t2;
        head_population->total++;
        //printf("\n\n ON %d population\n\n",i);
    }
    search_best_gene();
}
void printstats1(popptr t2){
    printstats(t2->pmstats,t2->vmstats,t2->cnstats,t2->pmvm,t2->vmcn,t2->id_counter,data1.cont_syn200);
}
void search_best_gene(){
    popptr t1=head_population->next;
    float best_utlz=0;
    popptr best;
    //remove->next here.
    while(t1->next!=NULL){
        printf("\nUtilization:%f",t1->utlz);
        if(best_utlz<t1->utlz){
            best_utlz=t1->utlz;
            best=t1;
        }
        t1=t1->next;
    }
    printf("\nbest ultilized is:%f",best_utlz);
};
void fitness(popptr pop){
    float energy[pop->id_counter[0]];
    float tot_utlz=0;
    for(int i=0;i<pop->id_counter[0];i++){
            energy[i]=0;
            int pid=i;
            for(int j=0;j<pop->id_counter[4];j++){
                if((pop->pmvm[j][0]==pid)&&(pop->pmvm[j][1]!=-1)){
                    int vid=pop->pmvm[j][1];
                    int vm_oid=pop->vmstats[j][2];
                    float vm_overhead=data1.vm_syn[vm_oid][0]*.1;
                    float cont_cpu=0;
                    for(int k=0;k<pop->id_counter[3];k++){
                        if((pop->vmcn[k][0]==vid)&&(pop->vmcn[k][1]!=-1)){
                            int cid=pop->vmcn[k][1];
                            int c_oid=data1.test_containers[pop->cnstats[cid]];
                            cont_cpu+=data1.cont_syn200[c_oid][0];
                        }
                    }
                    energy[i]+=vm_overhead+cont_cpu;
                }
            }
            energy[i]=energy[i]/data1.pm[0];
            pop->pmstats[i][3]=energy[i];
            printf("\nutilization of pm %d is %f",i,energy[i]);
            tot_utlz+=energy[i];
    }
    pop->utlz=tot_utlz/pop->id_counter[0];
    printf("\n\tTotal utilization is %f\n",pop->utlz);
}

popptr binary_selection(){
    popptr selected1=(popptr)malloc(sizeof(struct population));
    popptr selected2=(popptr)malloc(sizeof(struct population));

    int tt1=rand()%(head_population->total-1);
    int tt2=rand()%(head_population->total-1);
    while(tt2==tt1){
        tt2=rand()%head_population->total;
    }
    //printf("\n%d,%d\n",tt1,tt2);
    int i=0;
    float tf1,tf2;
    popptr t1=head_population->next;
    while(t1->next!=NULL){
        if(i==tt1){
            selected1=t1;
            tf1=t1->utlz;
        }
        if(i==tt2){
            selected2=t1;
            tf2=t1->utlz;
        }
        i++;
        t1=t1->next;
    }
    //printf("\ni:%d\n",i);
    //printf("\nUtlz:%f\n",selected1->utlz);
    //printf("\nUtlz:%f\n",selected2->utlz);
    if(tf2>tf1){
        popptr tempselected=selected2;
        selected2=selected1;
        selected1=tempselected;
    }

    int decider=rand()%10;

    if(decider>2){
        //printf("\nUtlz:%f\n",selected1->utlz);
        return(selected1);
    }
    else{
        //printf("\nUtlz:%f\n",selected2->utlz);
        return(selected2);
    }
}
void generate_random_containers(){
    int i=0;
    int check_bucket[200];
    while(i<200){
        check_bucket[i]=0;
        i++;
    }
    i=0;
    while(i<100){
        int ran=(rand()%200);
        while(check_bucket[ran]!=0){
            ran=(rand()%200);
        }
        data1.test_containers[i]=ran;
        check_bucket[ran]=1;
        i++;
    }
   }
void generate_random_containers2(int test_containers[100]){
    int i=0;
    int check_bucket[100];
    while(i<100){
        check_bucket[i]=0;
        i++;
    }
    i=0;
    while(i<100){
        int ran=(rand()%100);
        while(check_bucket[ran]!=0){
            ran=(rand()%100);
        }
        test_containers[i]=ran;
        check_bucket[ran]=1;
        //printf("\n%d added in test_containers\n",ran);
        i++;
    }
   }
void readcsv(){

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
        fscanf(fpr,"%f,%f\n",&data1.cont_syn200[i][0],&data1.cont_syn200[i][1]);
        i++;
    }
    //reading container file for real vms
    i=0;
    while(i<200){
        fscanf(fpr2,"%f,%f\n",&data1.cont_rea200[i][0],&data1.cont_rea200[i][1]);
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

void first_fit(
float pmstats[100][4],
float vmstats[100][4],
int cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float cont[200][2],
float vm[10][2],
float pm[3],
int NO_OF_TEST_CONTAINERS,
int test_containers[]
){

    for(int i=0;i<NO_OF_TEST_CONTAINERS;i++){
        int current_container=data1.test_containers[test_containers[i]];

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
            cnstats[id_counter[2]]=test_containers[i];
            int cid=id_counter[2];
            int vmcnid=id_counter[3];

            vmstats[selected_vid][0]-=cont[current_container][0];
            vmstats[selected_vid][1]-=cont[current_container][1];
            vmstats[selected_vid][3]++;

            vmcn[vmcnid][0]=selected_vid;
            vmcn[vmcnid][1]=cid;


            id_counter[2]++;
            id_counter[3]++;
            }
        else{
            //select random vm
            int selected_vm=(rand()%10);
            while(( cont[current_container][0] > vm[selected_vm][0] ) || ( cont[current_container][1] > vm[selected_vm][1] ) || ( vm[selected_vm][0] * 1.1 > pm[0] ) || (( vm[selected_vm][1] + 200 ) > pm[1] )){
                selected_vm=(rand()%10);
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

            cnstats[id_counter[2]]=test_containers[i];

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
float pmstats[100][4],
float vmstats[100][4],
int cnstats[100],
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
        int current_container=data1.test_containers[test_containers[i]];

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
            cnstats[id_counter[2]]=test_containers[i];
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

            cnstats[id_counter[2]]=test_containers[i];

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
float pmstats[100][4],
float vmstats[100][4],
int cnstats[100],
int pmvm[100][2],
int vmcn[100][2],
int id_counter[5],
float container[200][2]
){
    //for(int i=0;i<id_counter[3];)
    for(int i=0;i<id_counter[0];i++){
        printf("\n PM:%d Cpu left:%f Mem left:%f No. of VMs:%f VM:\n",i,pmstats[i][0],pmstats[i][1],pmstats[i][2]);
        for(int j=0;j<id_counter[4];j++){
            if((pmvm[j][0]==i)&&(pmvm[j][1]!=-1)){
                printf("\n\t VM:%d|%f Cpu left:%f Mem left:%f No. of Containers:%f Containers:\n",pmvm[j][1],vmstats[pmvm[j][1]][2],vmstats[pmvm[j][1]][0],vmstats[pmvm[j][1]][1],vmstats[pmvm[j][1]][3]);
                for(int k=0;k<id_counter[3];k++){
                    if((vmcn[k][0]==pmvm[j][1])&&(vmcn[k][1]!=-1)){
                        int cont=data1.test_containers[cnstats[vmcn[k][1]]];
                        printf("\n\t\t Container:%d||%d Cpu used:%f Mem used:%f.",vmcn[k][1],cnstats[vmcn[k][1]],container[cont][0],container[cont][1]);
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
