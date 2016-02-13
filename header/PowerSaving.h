#include<stdio.h> 
#include<time.h>
#include<stdlib.h>
#include"2Dquicksort4.h"
#include"quicksort2.h"
#include"3Dquicksort.h"

//#define SBP 8.00	//Power for BS in sleep mode
//#define ABP 68.73	//Basic power for BS in active mode

double * CoverageFinder(int *dBM, double *D, int Tbs, int Ntp){
	int i,j;
	double *coverage=malloc(Tbs*sizeof(double));
	ArrayInitialization(coverage,Tbs);
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			if(*(dBM+i*Ntp+j)==1)
				if(((*(D+(i+1)*(Ntp+1)+j+1))/20)>coverage[i])
					coverage[i]=(*(D+(i+1)*(Ntp+1)+j+1))/20;
					//printf("%g\n",coverage[i]);
	return coverage;
}

dBMprinter(int *dBM, int Ntp, int Tbs, int *DS, double *P){
               int i,j;
               for(i=0;i<Tbs;i++)
               printf("       BS%3d",i+1);
               printf("\n");
               
               for(i=0;i<Ntp;i++){
                                  printf("MS%4d%4d",i+1,*(dBM+i));
                                  for(j=1;j<Tbs;j++)
                                  printf("      %6d",*(dBM+j*Ntp+i));
                                  printf("\n");
                                  
                                  printf(" Power%10.5f",*(P+(Ntp+1)+i+1));
                                  for(j=1;j<Tbs;j++)
                                  printf("  %10.5f",*(P+(j+1)*(Ntp+1)+i+1));
                                  printf("\n");
                                  
                                  printf(" DS %6d",*(DS+i+1));
                                  for(j=1;j<Tbs;j++)
                                  printf("      %6d",*(DS+i+1));
                                  printf("\n");
                                  }
                                  }

temp_dBM_Initiator(int *ptr, int Ntp){	//initiate the array of temp_dBM for the new heuristic algorithm
	int i;
	for(i=0;i<Ntp;i++)
		*(ptr+i)=-1;
}

ArrayInitialization(double *ptr, int Tbs){
       int i;         
       for(i=0;i<Tbs;i++)
                          *(ptr+i)=0;

       
}
ArrayInitialization1(int *ptr, int Tbs){
       int i;         
       for(i=0;i<Tbs;i++)
                          *(ptr+i)=0;

       
}     
ArrayInitialization2(int *ptr, int Tbs, int Ntp){
       int i,j;         
       for(i=0;i<Tbs;i++)
       for(j=0;j<Ntp;j++)
       *(ptr+i*Ntp+j)=0;

       
}           
heuristic(int Tbs, double *P, int Ntp, int MP, int DSt, double BP, int *DS, int *Xbs, int *Ybs, int *Xtp, int *Ytp, double *ptrD, double *BR, int Modulation, double mABP, double mSBP){
                
       double *SortedPowerRequirement;  
       int *SortedPowerRequirementIndex;
       
       double UnsortedCumulativePower[Tbs];
       ArrayInitialization(UnsortedCumulativePower,Tbs);
       double SortedCumulativePower[Tbs];
       double TotalPower=0;
	double NOV=0.0;
       int *SortedCumulativePowerIndex;
       
       int CumulativeDS[Tbs];
       ArrayInitialization1(CumulativeDS,Tbs);
       int TotalDSs=0;
                            
       int dBM[Tbs][Ntp];
       ArrayInitialization2(dBM[0],Tbs,Ntp);
       int CumulativeServingMSs[Tbs];
       ArrayInitialization1(CumulativeServingMSs,Tbs);
       int TotalServedMSs=0;
       
	int BStatus[Tbs];
	ArrayInitialization1(BStatus,Tbs);
	
	double Bandwidth[Tbs];
	ArrayInitialization(Bandwidth,Tbs);
	double TotalBandwidth;

       int i,j,k,l;
       double MaxP=0.0;
       int IndexMaxP=-1;

       double *coverage;
           
       clock_t start_tick, end_tick;
       double elapsed;
       
       start_tick=clock();
	double *PIndex;
	int m,n;
       SortedPowerRequirement=SortedPowerRequirementAssignment(P,Tbs,Ntp);
       PIndex=q_sort_3D(P,Tbs,Ntp,Modulation);
       for(n=0;n<Ntp;n++){
		for(m=0;m<Tbs*Modulation;m++){
			printf("%d MS%d and BS%d use modulation%d transmitting power %lf. \n",n*Tbs*Modulation+m,n,(int)floor(*(PIndex+m*Ntp+n)),(int)(*(PIndex+m*Ntp+n)*10-floor(*(PIndex+m*Ntp+n))*10),*(P+((int)floor(*(PIndex+m*Ntp+n)))*Ntp*Modulation+n*Modulation+(int)(*(PIndex+m*Ntp+n)*10-floor(*(PIndex+m*Ntp+n))*10)));
			//printf("%d MS%d %lf transmitting power %lf. \n",n*Tbs*Modulation+m,n,(*(IndexBSMod+n*(Tbs*Modulation)+m)),*(BSMod+n*(Tbs*Modulation)+m));
		}
		printf("\n");
	}
                         
       SortedPowerRequirementIndex=q_sort_2D(SortedPowerRequirement,0,Tbs-1,Ntp,Tbs);
       
      
       /*for(k=0;k<Ntp;k++){
       for(j=0;j<Tbs;j++)
       printf("Sorted Power Requirement from BS %d to MS %d: %lf \n",*(SortedPowerRequirementIndex+j*Ntp+k)+1,k+1,*(SortedPowerRequirement+j*Ntp+k));
       printf("\n");
       }*/
       
       //new heuristic algorithm goes from here
       int temp_dBM[Ntp];
       temp_dBM_Initiator(temp_dBM,Ntp);
	double MinP=1000.0;
	int IndexMinP=-1;	

       for(i=0;i<Ntp;i++){	//each MS is assigned to the 1st BS choice
	*(temp_dBM+i)=*(SortedPowerRequirementIndex+i);
	CumulativeServingMSs[*(SortedPowerRequirementIndex+i)]++;
	TotalServedMSs++;
	*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i))+=*(SortedPowerRequirement+i);
	*(CumulativeDS+*(SortedPowerRequirementIndex+i))+=*(DS+i+1);
	BStatus[*(SortedPowerRequirementIndex+i)]=1;
}
int temp=0;
for(i=0;i<Tbs;i++){
	temp+=CumulativeServingMSs[i];
	//printf("%d\n", i);
}
if(TotalServedMSs!=temp){
printf("0. Inconsistent numbers! %d %d %d %d %d %d",TotalServedMSs,CumulativeServingMSs[0],CumulativeServingMSs[1],CumulativeServingMSs[2],CumulativeServingMSs[3],temp);
exit(1);}

for(j=0;j<Tbs;j++){
	while((*(UnsortedCumulativePower+j)>MP)||(*(CumulativeDS+j)>DSt)){
		for(i=0;i<Ntp;i++){
			if((*(temp_dBM+i)==j)&&(*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1)<MinP)){
				MinP=*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1);
				IndexMinP=i;
			}
		}
		if(IndexMinP!=-1){
			*(UnsortedCumulativePower+*(temp_dBM+IndexMinP))-=*(P+(*(temp_dBM+IndexMinP)+1)*(Ntp+1)+IndexMinP+1);
			*(CumulativeDS+*(temp_dBM+IndexMinP))-=*(DS+IndexMinP+1);
			CumulativeServingMSs[*(temp_dBM+IndexMinP)]--;
			TotalServedMSs--;
			*(temp_dBM+IndexMinP)=-1;
		}
	BStatus[j]=2;
	IndexMinP=-1;
	MinP=1000.0;
	}	
}

if(TotalServedMSs<ceil(Ntp*(1-BP))){
	for(i=0;i<Ntp;i++){
		if(*(temp_dBM+i)==-1){
			*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i+Ntp))+=*(SortedPowerRequirement+i+Ntp);
			*(CumulativeDS+*(SortedPowerRequirementIndex+i+Ntp))+=*(DS+i+1);
			*(temp_dBM+i)=*(SortedPowerRequirementIndex+i+Ntp);
			CumulativeServingMSs[*(SortedPowerRequirementIndex+i+Ntp)]++;
			TotalServedMSs++;
			BStatus[*(SortedPowerRequirementIndex+i+Ntp)]=1;
		}
	}
}

for(j=0;j<Tbs;j++){
	while((*(UnsortedCumulativePower+j)>MP)||(*(CumulativeDS+j)>DSt)){
		for(i=0;i<Ntp;i++){
			if((*(temp_dBM+i)==j)&&(*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1)>MaxP)){
				MaxP=*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1);
				IndexMaxP=i;
			}
		}
		if(IndexMaxP!=-1){
			*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))-=*(P+(*(temp_dBM+IndexMaxP)+1)*(Ntp+1)+IndexMaxP+1);
			*(CumulativeDS+*(temp_dBM+IndexMaxP))-=*(DS+IndexMaxP+1);
			CumulativeServingMSs[*(temp_dBM+IndexMaxP)]--;
			TotalServedMSs--;
			*(temp_dBM+IndexMaxP)=-1;
		}
	BStatus[j]=2;
	IndexMaxP=-1;
	MaxP=0.0;
	}	
}
/*printf("TotalServedMSs %d\n",TotalServedMSs);
for(i=0;i<Tbs;i++)
	printf("CumulativeServingMSs %d\n",CumulativeServingMSs[i]);*/
/*for(i=0;i<Tbs;i++)
	if(BStatus[i]==0)
		*(UnsortedCumulativePower+i)+=SBP;
	else
		*(UnsortedCumulativePower+i)+=ABP;*/
for(i=0;i<Tbs;i++)
       *(SortedCumulativePower+i)=*(UnsortedCumulativePower+i);
       
//        printf("\n");
       SortedCumulativePowerIndex=DOUBLE_q_sort(SortedCumulativePower,0,Tbs-1,Tbs);
/*       for(i=0;i<Tbs;i++)
       printf("Unsorted Cumulative Power of BS %d: %lf \n",i+1,*(UnsortedCumulativePower+i));
       printf("\n");
       for(i=0;i<Tbs;i++)
       printf("%d. Sorted Cumulative Power of BS %d: %lf \n",i+1,*(SortedCumulativePowerIndex+i)+1,*(SortedCumulativePower+i));
       printf("\n");*/
       /*for(i=0;i<Tbs;i++)
       printf("%d ",*(SortedCumulativePowerIndex+i)+1);
       printf("\nTotal Served Users:%d\n",TotalServedMSs);
	for(i=0;i<Tbs;i++)
		printf("BS%d is in status %d.\n",i+1,BStatus[i]);*/
/*for(l=0;l<Tbs;l++)
	if(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt){
		printf("Program Stops!\n");
		exit(1);
}*/

double difference;
int RemovedUsers=0;
	for(i=Tbs-1;i>-1;i--){
		//while(difference>0&&difference<ABP-SBP)
			if(BStatus[*(SortedCumulativePowerIndex+i)]==1){
				difference=0.0;
				RemovedUsers=0;
				double *AdditionPower=malloc(Tbs*sizeof(double));
				int *AdditionDS=malloc(Tbs*sizeof(int));
				//double AdditionPower[Tbs];
				//int AdditionDS[Tbs];
				ArrayInitialization(AdditionPower,Tbs);
				ArrayInitialization1(AdditionDS,Tbs);
				for(j=0;j<Ntp;j++){
					if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i))
						for(k=0;k<Tbs;k++)
							if((BStatus[*(SortedPowerRequirementIndex+k*Ntp+j)]==1)&&(*(SortedPowerRequirementIndex+k*Ntp+j)!=*(SortedCumulativePowerIndex+i))&&*(AdditionDS+*(SortedPowerRequirementIndex+k*Ntp+j))+(*(DS+j+1))+(*(CumulativeDS+*(SortedPowerRequirementIndex+k*Ntp+j)))<=DSt&&*(AdditionPower+*(SortedPowerRequirementIndex+k*Ntp+j))+(*(P+(*(SortedPowerRequirementIndex+k*Ntp+j)+1)*(Ntp+1)+j+1))+(*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+k*Ntp+j)))<=MP){
								difference+=*(P+(*(SortedPowerRequirementIndex+k*Ntp+j)+1)*(Ntp+1)+j+1)-*(P+(*(SortedCumulativePowerIndex+i)+1)*(Ntp+1)+j+1);
								
								*(AdditionPower+*(SortedPowerRequirementIndex+k*Ntp+j))+=*(P+(*(SortedPowerRequirementIndex+k*Ntp+j)+1)*(Ntp+1)+j+1);
								*(AdditionDS+*(SortedPowerRequirementIndex+k*Ntp+j))+=*(DS+j+1);
								//if(AddtionPower[*(SortedPowerRequirementIndex+k*Ntp+j)]+)
								RemovedUsers++;
								k=Tbs;
							}
				}
				free(AdditionPower);
				free(AdditionDS);
				if(difference>0&&difference<ABP-SBP&&RemovedUsers==*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i))){
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i))
							for(k=0;k<Tbs;k++)
								if(BStatus[*(SortedPowerRequirementIndex+k*Ntp+j)]==1&&*(SortedPowerRequirementIndex+k*Ntp+j)!=*(SortedCumulativePowerIndex+i)&&(*(CumulativeDS+*(SortedPowerRequirementIndex+k*Ntp+j)))+(*(DS+j+1))<=DSt&&((*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+k*Ntp+j)))+(*(P+(*(SortedPowerRequirementIndex+k*Ntp+j)+1)*(Ntp+1)+j+1)))<=MP){
									*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+k*Ntp+j))+=*(P+(*(SortedPowerRequirementIndex+k*Ntp+j)+1)*(Ntp+1)+j+1);
// 									printf("before %lf\n",*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
									*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))-=*(P+(*(SortedCumulativePowerIndex+i)+1)*(Ntp+1)+j+1);
// 									printf("after %lf, substracted by %lf\n",*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)),*(P+(*(SortedCumulativePowerIndex+i)+1)*(Ntp+1)+j+1));
									*(CumulativeDS+*(SortedPowerRequirementIndex+k*Ntp+j))+=*(DS+j+1);
									*(CumulativeDS+*(SortedCumulativePowerIndex+i))-=*(DS+j+1);
									*(CumulativeServingMSs+*(SortedPowerRequirementIndex+k*Ntp+j))+=1;
									*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i))-=1;
									*(temp_dBM+j)=*(SortedPowerRequirementIndex+k*Ntp+j);
									if((*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+k*Ntp+j)))==MP||(*(CumulativeDS+*(SortedPowerRequirementIndex+k*Ntp+j))==DSt))
											BStatus[*(SortedPowerRequirementIndex+k*Ntp+j)]=2;
									k=Tbs;
									
								}
					if((*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i)))==0){
						BStatus[*(SortedCumulativePowerIndex+i)]=0;
						//*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))-=ABP+SBP;
					}
					else{
						printf("4. Program Stops! BS%d serves %d users instead of 0.\n",*(SortedCumulativePowerIndex+i)+1,*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i)));
						exit(1);
					}
						
					for(l=0;l<Tbs;l++)
						if(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt){
							printf("Program Stops! BS%d transmits power at %lf level.\n",l+1,*(UnsortedCumulativePower+l));
							exit(1);
							//BStatus[l]=2;
						}
					for(l=0;l<Tbs;l++)
						*(SortedCumulativePower+l)=*(UnsortedCumulativePower+l);
					SortedCumulativePowerIndex=DOUBLE_q_sort(SortedCumulativePower,0,Tbs-1,Tbs);

					
					i=Tbs;
				}
				
			}
				
	}
/*for(l=0;l<Tbs;l++)
	if(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt){
		printf("Program Stops!\n");
		exit(1);
}*/

double NewTransmissionPower;
double MaxSavablePower;
double SavablePower;
int NewBS;
int UserCutOff;
for(i=Tbs-1;i>-1;i--){
	NewBS=-1;
	MaxSavablePower=0.0;
	UserCutOff=-1;
	if(BStatus[*(SortedCumulativePowerIndex+i)]==1||BStatus[*(SortedCumulativePowerIndex+i)]==2){
		//printf("Old BS:%d\n",*(SortedCumulativePowerIndex+i));
		for(k=0;k<Tbs;k++){
			if(*(SortedCumulativePowerIndex+i)!=*(SortedCumulativePowerIndex+k)){
				NewTransmissionPower=0.0;
				SavablePower=0.0;
				//printf("Old BS:%d, New BS:%d\n",*(SortedCumulativePowerIndex+i)+1,*(SortedCumulativePowerIndex+k)+1);
				if((BStatus[*(SortedCumulativePowerIndex+k)]==1)||(BStatus[*(SortedCumulativePowerIndex+k)]==0))
					//if(((*(CumulativeDS+*(SortedCumulativePowerIndex+k)))+(*(CumulativeDS+*(SortedCumulativePowerIndex+i))))<=DSt)
						for(j=0;j<Ntp;j++)
							if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i))
								NewTransmissionPower+=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
				//printf("To New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",*(SortedCumulativePowerIndex+k)+1,NewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
				NewTransmissionPower+=(*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k)))+ABP+SBP;
				if(BStatus[*(SortedCumulativePowerIndex+k)]==0){	
					if(NewTransmissionPower-ABP-SBP<=MP){
						SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+SBP-NewTransmissionPower;
						if(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))<=DSt){
							MaxSavablePower=SavablePower;
							NewBS=*(SortedCumulativePowerIndex+k);
							UserCutOff=0;
						}
					}
					else
					if(NewTransmissionPower-ABP-SBP>MP){
						int CombinedSize=*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))+*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i));
						double *CombinedBSs=calloc(CombinedSize,sizeof(double));
						int *CombinedDSs=calloc(CombinedSize,sizeof(int));
						int RemovableUsers=TotalServedMSs-ceil(Ntp*(1-BP));
						int *Index;
						int XDS=0;
						double Xpower=0.0;
						int m=0;
						if(RemovableUsers>0){
							for(j=0;j<Ntp;j++)
								if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){
									CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
									CombinedDSs[m]=*(DS+j+1);
									m++;
								}
							for(j=0;j<Ntp;j++)
								if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+k)){
									CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
									CombinedDSs[m]=*(DS+j+1);
									m++;
								}
							Index=DOUBLE_q_sort(CombinedBSs,0,CombinedSize-1,CombinedSize);
							j=CombinedSize-1;
							if(CombinedSize>RemovableUsers)
								while(j>CombinedSize-RemovableUsers-1){
									Xpower+=CombinedBSs[j];
									XDS+=CombinedDSs[*(Index+j)];
									SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+SBP-NewTransmissionPower+Xpower;
									if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
										MaxSavablePower=SavablePower;
										NewBS=*(SortedCumulativePowerIndex+k);
										UserCutOff=1;
										j=-1;
	
									}
									else
										j--;
								}
							else if(CombinedSize<=RemovableUsers){
								while(j>-1){
									Xpower+=CombinedBSs[j];
									XDS+=CombinedDSs[*(Index+j)];
									SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+SBP-NewTransmissionPower+Xpower;
									if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
										MaxSavablePower=SavablePower;
										NewBS=*(SortedCumulativePowerIndex+k);
										UserCutOff=1;
										j=-1;
	
									}
									else
										j--;
								}
							}
							
							
						}
						free(CombinedBSs);
						free(CombinedDSs);
					}
				}	
				else 
				if(BStatus[*(SortedCumulativePowerIndex+k)]==1){	
					if(NewTransmissionPower-ABP-SBP<=MP){
						SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower;
						if(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))<=DSt){
							MaxSavablePower=SavablePower;
							NewBS=*(SortedCumulativePowerIndex+k);
							UserCutOff=0;
						}
						else if(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))>DSt){
							int CombinedSize=*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))+*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i));
							double *CombinedBSs=calloc(CombinedSize,sizeof(double));
							double *CombinedDSs=calloc(CombinedSize,sizeof(double));
							int RemovableUsers=TotalServedMSs-ceil(Ntp*(1-BP));
							int *Index;
							int XDS=0;
							double Xpower=0.0;
							int m=0;
							if(RemovableUsers>0){
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){
										CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
										CombinedDSs[m]=*(DS+j+1);
										m++;
									}
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+k)){
										CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
										CombinedDSs[m]=*(DS+j+1);
										m++;
									}
								Index=DOUBLE_q_sort(CombinedDSs,0,CombinedSize-1,CombinedSize);
								j=CombinedSize-1;
								if(CombinedSize>RemovableUsers)
									while(j>CombinedSize-RemovableUsers-1){
										XDS+=CombinedDSs[j];
										Xpower+=CombinedBSs[*(Index+j)];
										SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
										if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
											MaxSavablePower=SavablePower;
											NewBS=*(SortedCumulativePowerIndex+k);
											UserCutOff=2;
											j=-1;
										}
										else
											j--;
									}
								else if(CombinedSize<=RemovableUsers)
									while(j>-1){
										XDS+=CombinedDSs[j];
										Xpower+=CombinedBSs[*(Index+j)];
										SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
										if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
											MaxSavablePower=SavablePower;
											NewBS=*(SortedCumulativePowerIndex+k);
											UserCutOff=2;
											j=-1;
										}
										else
											j--;
									}
								
							free(CombinedBSs);
							free(CombinedDSs);
							}
						}	
					}	
					else 
					if(NewTransmissionPower-ABP-SBP>MP){
						if(*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))<=DSt){
							int CombinedSize=*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))+*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i));
							double *CombinedBSs=calloc(CombinedSize,sizeof(double));
							int *CombinedDSs=calloc(CombinedSize,sizeof(int));
							int RemovableUsers=TotalServedMSs-ceil(Ntp*(1-BP));
							int *Index;
							int XDS=0;
							double Xpower=0.0;
							int m=0;
							if(RemovableUsers>0){
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){
										CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
										CombinedDSs[m]=*(DS+j+1);
										m++;
									}
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+k)){
										CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
										CombinedDSs[m]=*(DS+j+1);
										m++;
									}
								Index=DOUBLE_q_sort(CombinedBSs,0,CombinedSize-1,CombinedSize);
								j=CombinedSize-1;
								if(CombinedSize>RemovableUsers)
									while(j>CombinedSize-RemovableUsers-1){
										Xpower+=CombinedBSs[j];
										XDS+=CombinedDSs[*(Index+j)];
										SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
										if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
											MaxSavablePower=SavablePower;
											NewBS=*(SortedCumulativePowerIndex+k);
											UserCutOff=1;
											j=-1;
		
										}
										else
											j--;
									}
								else if(CombinedSize<=RemovableUsers)
									while(j>-1){
										Xpower+=CombinedBSs[j];
										XDS+=CombinedDSs[*(Index+j)];
										SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
										if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
											MaxSavablePower=SavablePower;
											NewBS=*(SortedCumulativePowerIndex+k);
											UserCutOff=1;
											j=-1;
		
										}
										else
											j--;
									}
								
								
							}
							free(CombinedBSs);
							free(CombinedDSs);
						}
						else 
						if(*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))>DSt){
							double DSRatio=(*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-DSt)/DSt;
							double PRatio=(NewTransmissionPower-ABP-SBP-MP)/MP;
							if(DSRatio>PRatio){
								int CombinedSize=*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))+*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i));
								double *CombinedBSs=calloc(CombinedSize,sizeof(double));
								double *CombinedDSs=calloc(CombinedSize,sizeof(double));
								int RemovableUsers=TotalServedMSs-ceil(Ntp*(1-BP));
								int *Index;
								int XDS=0;
								double Xpower=0.0;
								int m=0;
								if(RemovableUsers>0){
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){
											CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
											CombinedDSs[m]=*(DS+j+1);
											m++;
										}
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+k)){
											CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
											CombinedDSs[m]=*(DS+j+1);
											m++;
										}
									Index=DOUBLE_q_sort(CombinedDSs,0,CombinedSize-1,CombinedSize);
									j=CombinedSize-1;
									if(CombinedSize>RemovableUsers)
										while(j>CombinedSize-RemovableUsers-1){
											XDS+=CombinedDSs[j];
											Xpower+=CombinedBSs[*(Index+j)];
											SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
											if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
												MaxSavablePower=SavablePower;
												NewBS=*(SortedCumulativePowerIndex+k);
												UserCutOff=2;
												j=-1;
											}
											else
												j--;
										}
									else if(CombinedSize<=RemovableUsers)
										while(j>-1){
											XDS+=CombinedDSs[j];
											Xpower+=CombinedBSs[*(Index+j)];
											SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
											if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
												MaxSavablePower=SavablePower;
												NewBS=*(SortedCumulativePowerIndex+k);
												UserCutOff=2;
												j=-1;
											}
											else
												j--;
										}
									
								free(CombinedBSs);
								free(CombinedDSs);
								}
							}
								
							else if(DSRatio<=PRatio){
								int CombinedSize=*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))+*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i));
								double *CombinedBSs=calloc(CombinedSize,sizeof(double));
								int *CombinedDSs=calloc(CombinedSize,sizeof(int));
								int RemovableUsers=TotalServedMSs-ceil(Ntp*(1-BP));
								int *Index;
								int XDS=0;
								double Xpower=0.0;
								int m=0;
								if(RemovableUsers>0){
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){
											CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
											CombinedDSs[m]=*(DS+j+1);
											m++;
										}
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+k)){
											CombinedBSs[m]=*(P+(*(SortedCumulativePowerIndex+k)+1)*(Ntp+1)+j+1);
											CombinedDSs[m]=*(DS+j+1);
											m++;
										}
									Index=DOUBLE_q_sort(CombinedBSs,0,CombinedSize-1,CombinedSize);
									j=CombinedSize-1;
									if(CombinedSize>RemovableUsers)
										while(j>CombinedSize-RemovableUsers-1){
											Xpower+=CombinedBSs[j];
											XDS+=CombinedDSs[*(Index+j)];
											SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
											if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
												MaxSavablePower=SavablePower;
												NewBS=*(SortedCumulativePowerIndex+k);
												UserCutOff=1;
												j=-1;
			
											}
											else
												j--;
										}
									else if(CombinedSize<=RemovableUsers)
										while(j>-1){
											Xpower+=CombinedBSs[j];
											XDS+=CombinedDSs[*(Index+j)];
											SavablePower=*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))+ABP+*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))+ABP-NewTransmissionPower+Xpower;
											if(NewTransmissionPower-ABP-SBP-Xpower<=MP&&(SavablePower>MaxSavablePower&&*(CumulativeDS+*(SortedCumulativePowerIndex+i))+*(CumulativeDS+*(SortedCumulativePowerIndex+k))-XDS<=DSt)){
												MaxSavablePower=SavablePower;
												NewBS=*(SortedCumulativePowerIndex+k);
												UserCutOff=1;
												j=-1;
			
											}
											else
												j--;
										}
									
									
								}
								free(CombinedBSs);
								free(CombinedDSs);
							}
						}
					}	
				}
			}		
			//printf("New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",NewBS+1,MinNewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
			
		}	
	}
// 	for(l=0;l<Tbs;l++)
// 		printf("BS%d is in status %d.\n",l+1,BStatus[l]);
	//printf("New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",NewBS+1,MinNewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
	if(NewBS!=-1&&UserCutOff==0){
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){							
				*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
				*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))-=*(P+(*(SortedCumulativePowerIndex+i)+1)*(Ntp+1)+j+1);
				*(CumulativeDS+NewBS)+=*(DS+j+1);
				*(CumulativeDS+*(SortedCumulativePowerIndex+i))-=*(DS+j+1);
				*(CumulativeServingMSs+NewBS)+=1;
				*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i))-=1;
				*(temp_dBM+j)=NewBS;
				
			}
		//printf("BS%d(%d) Swith to New BS:%d(%d)\n",*(SortedCumulativePowerIndex+i)+1,BStatus[*(SortedCumulativePowerIndex+i)],NewBS+1,BStatus[NewBS]);
		BStatus[*(SortedCumulativePowerIndex+i)]=0;
		if((*(UnsortedCumulativePower+NewBS))==MP||(*(CumulativeDS+NewBS)==DSt))
			BStatus[NewBS]=2;
		else
			BStatus[NewBS]=1;
		//printf("BS%d(%d) Swith to New BS:%d(%d)\n",*(SortedCumulativePowerIndex+i)+1,BStatus[*(SortedCumulativePowerIndex+i)],NewBS+1,BStatus[NewBS]);
		//for(l=0;l<Tbs;l++)
		//	printf("BS%d is in status %d.\n",l+1,BStatus[l]);
		i=Tbs;
		for(l=0;l<Tbs;l++)
			if(BStatus[l]!=0&&(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt)){
				printf("0. Program Stops! BS%d transmits power at %lf level.\n",l+1,*(UnsortedCumulativePower+l));
				exit(1);
			}
			else
			if(BStatus[l]==0&&(*(UnsortedCumulativePower+l)>0.00000000001||*(CumulativeDS+l)>0.000000000001)){
				printf("1. Program Stops! BS%d transmits power at %lf level.\n",l+1,*(UnsortedCumulativePower+l));
				exit(1);
			}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativePower+l)=*(UnsortedCumulativePower+l);
		SortedCumulativePowerIndex=DOUBLE_q_sort(SortedCumulativePower,0,Tbs-1,Tbs);
	}
	else 
	if(NewBS!=-1&&UserCutOff==1){
		
		if(MaxSavablePower>0){
// 		if(MinNewTransmissionPower+*(UnsortedCumulativePower+NewBS)-SupressPower<=MP){
			for(j=0;j<Ntp;j++)
				if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){
					*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
					*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))-=*(P+(*(SortedCumulativePowerIndex+i)+1)*(Ntp+1)+j+1);
					*(CumulativeDS+NewBS)+=*(DS+j+1);
					*(CumulativeDS+*(SortedCumulativePowerIndex+i))-=*(DS+j+1);
					*(CumulativeServingMSs+NewBS)+=1;
					*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i))-=1;
					*(temp_dBM+j)=NewBS;
				}
			if(*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i))==0)
				BStatus[*(SortedCumulativePowerIndex+i)]=0;
			else{
				printf("Program Stops! BS%d transmits power at %lf level.\n",*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
				exit(1);
			}
			IndexMaxP=-1;
			MaxP=0.0;
			while(*(UnsortedCumulativePower+NewBS)>MP||*(CumulativeDS+NewBS)>DSt){
				for(j=0;j<Ntp;j++){
					if(*(temp_dBM+j)==NewBS&&*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1)>MaxP){
						MaxP=*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1);
						IndexMaxP=j;
					}
				}
				if(IndexMaxP!=-1){
					*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))-=*(P+(*(temp_dBM+IndexMaxP)+1)*(Ntp+1)+IndexMaxP+1);
					*(CumulativeDS+*(temp_dBM+IndexMaxP))-=*(DS+IndexMaxP+1);
					CumulativeServingMSs[*(temp_dBM+IndexMaxP)]--;
					TotalServedMSs--;
					
					if(*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))==MP||(*(CumulativeDS+*(temp_dBM+IndexMaxP))==DSt))
						BStatus[*(temp_dBM+IndexMaxP)]=2;
					else
						BStatus[*(temp_dBM+IndexMaxP)]=1;
					*(temp_dBM+IndexMaxP)=-1;
				}
				IndexMaxP=-1;
				MaxP=0;
			}	
			
		i=Tbs;
		for(l=0;l<Tbs;l++)
			if(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt){
				printf("Program Stops! BS%d transmits power at %lf level.\n",l+1,*(UnsortedCumulativePower+l));
				exit(1);
			}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativePower+l)=*(UnsortedCumulativePower+l);
		SortedCumulativePowerIndex=DOUBLE_q_sort(SortedCumulativePower,0,Tbs-1,Tbs);
		}
		UserCutOff=0;
		
	}
	else 
	if(NewBS!=-1&&UserCutOff==2){
		if(MaxSavablePower>0){
			for(j=0;j<Ntp;j++)
				if(*(temp_dBM+j)==*(SortedCumulativePowerIndex+i)){
					*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
					*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i))-=*(P+(*(SortedCumulativePowerIndex+i)+1)*(Ntp+1)+j+1);
					*(CumulativeDS+NewBS)+=*(DS+j+1);
					*(CumulativeDS+*(SortedCumulativePowerIndex+i))-=*(DS+j+1);
					*(CumulativeServingMSs+NewBS)+=1;
					*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i))-=1;
					*(temp_dBM+j)=NewBS;
				}
			if(*(CumulativeServingMSs+*(SortedCumulativePowerIndex+i))==0)
				BStatus[*(SortedCumulativePowerIndex+i)]=0;
			else{
				printf("Program Stops! BS%d transmits power at %lf level.\n",*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
				exit(1);
			}
			int IndexMaxDS=-1;
			int MaxDS=0;
			while(*(UnsortedCumulativePower+NewBS)>MP||*(CumulativeDS+NewBS)>DSt){
				for(j=0;j<Ntp;j++){
					if(*(temp_dBM+j)==NewBS&&*(DS+j+1)>MaxDS){
						MaxDS=*(DS+j+1);
						IndexMaxDS=j;
					}
				}
				if(IndexMaxDS!=-1){
					*(UnsortedCumulativePower+*(temp_dBM+IndexMaxDS))-=*(P+(*(temp_dBM+IndexMaxDS)+1)*(Ntp+1)+IndexMaxDS+1);
					*(CumulativeDS+*(temp_dBM+IndexMaxDS))-=*(DS+IndexMaxDS+1);
					CumulativeServingMSs[*(temp_dBM+IndexMaxDS)]--;
					TotalServedMSs--;
					
					if(*(UnsortedCumulativePower+*(temp_dBM+IndexMaxDS))==MP||(*(CumulativeDS+*(temp_dBM+IndexMaxDS))==DSt))
						BStatus[*(temp_dBM+IndexMaxDS)]=2;
					else
						BStatus[*(temp_dBM+IndexMaxDS)]=1;
					*(temp_dBM+IndexMaxDS)=-1;
				}
				IndexMaxDS=-1;
				MaxDS=0;
			}	
			
		i=Tbs;
		for(l=0;l<Tbs;l++)
			if(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt){
				printf("Program Stops! BS%d transmits power at %lf level.\n",l+1,*(UnsortedCumulativePower+l));
				exit(1);
			}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativePower+l)=*(UnsortedCumulativePower+l);
		SortedCumulativePowerIndex=DOUBLE_q_sort(SortedCumulativePower,0,Tbs-1,Tbs);
		}

		
	}		
		
		
	
}		
				
IndexMaxP=-1;
MaxP=0;
while(TotalServedMSs>ceil(Ntp*(1-BP))){
	//for(i=0;i<Ntp;i++)		printf("MS%d, index: %d, power: %lf\n ",i+1,*(temp_dBM+i)+1,*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1));
	for(i=0;i<Ntp;i++){
		if(*(temp_dBM+i)!=-1&&*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1)>MaxP){
			MaxP=*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1);
			IndexMaxP=i;
		}
	}
	if(IndexMaxP!=-1){
		*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))-=*(P+(*(temp_dBM+IndexMaxP)+1)*(Ntp+1)+IndexMaxP+1);
		*(CumulativeDS+*(temp_dBM+IndexMaxP))-=*(DS+IndexMaxP+1);
		CumulativeServingMSs[*(temp_dBM+IndexMaxP)]--;
		TotalServedMSs--;
//printf("BS %d SERVES # %d\n",*(temp_dBM+IndexMaxP),CumulativeServingMSs[*(temp_dBM+IndexMaxP)]);
		
		if(CumulativeServingMSs[*(temp_dBM+IndexMaxP)]==0)
			BStatus[*(temp_dBM+IndexMaxP)]=0;
		else
			BStatus[*(temp_dBM+IndexMaxP)]=1;
		*(temp_dBM+IndexMaxP)=-1;
	}
	IndexMaxP=-1;
	MaxP=0;
	//printf("REMOVE INDEX%d, %lf",IndexMaxP,*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP)));
}	

temp=0;
/*for(i=0;i<Tbs;i++)
	printf("BS%d is in status %d.\n",i+1,BStatus[i]);*/
	//temp+=CumulativeServingMSs[i];
/*if(TotalServedMSs!=temp){
printf("4. Inconsistent numbers!%d %d",temp,TotalServedMSs);
exit(1);}*/

for(j=0;j<Tbs;j++){
	if((*(UnsortedCumulativePower+j)>MP)||(*(CumulativeDS+j)>DSt)){
		printf("Constraint violation! Program Stops!!\n");
		exit(1);
	}
}

	for(i=0;i<Ntp;i++)
		if(*(temp_dBM+i)!=-1){
			dBM[*(temp_dBM+i)][i]=1;
			Bandwidth[*(temp_dBM+i)]+=*(BR+i+1);
			TotalBandwidth+=*(BR+i+1);
			
		}
temp=0;
for(i=0;i<Tbs;i++)
	temp+=CumulativeServingMSs[i];
if(TotalServedMSs!=temp){
printf("6. Inconsistent numbers!%d %d",temp,TotalServedMSs);
exit(1);}
       /*for(i=0;i<Ntp;i++){
                          if(((*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i))+*(SortedPowerRequirement+i))<=MP)&&
                             (*(CumulativeDS+*(SortedPowerRequirementIndex+i))+*(DS+i+1))<=DSt){	//DSt*200 for 1 second
                              //printf("Before assignment, MS %d needs DS %d, BS cumulative DSs %d: %d\n",i,*(DS+i+1),*(SortedPowerRequirementIndex+i)+1,*(CumulativeDS+*(SortedPowerRequirementIndex+i)));
                              *(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i))+=*(SortedPowerRequirement+i);
                              *(CumulativeDS+*(SortedPowerRequirementIndex+i))+=*(DS+i+1);
                              dBM[*(SortedPowerRequirementIndex+i)][i]=1;
                              CumulativeServingMSs[*(SortedPowerRequirementIndex+i)]++;
                              TotalServedMSs++;
                              //printf("After assignment, MS %d needs DS %d, BS cumulative DSs %d: %d\n",i,*(DS+i+1),*(SortedPowerRequirementIndex+i)+1,*(CumulativeDS+*(SortedPowerRequirementIndex+i)));
                          }
                          else{
			       for(k=0;k<=Ntp;k++)
                               if((*(dBM+*(SortedPowerRequirementIndex+i)*Ntp+k)==1)&&(*(P+*(SortedPowerRequirementIndex+i)*Ntp+k)>MaxP)){
                               MaxP=*(P+*(SortedPowerRequirementIndex+i)*Ntp+k);
                               IndexMaxP=k;
                               }
                               if((*(P+*(SortedPowerRequirementIndex+i)*Ntp+IndexMaxP)>*(SortedPowerRequirement+i))&&((*(CumulativeDS+*(SortedPowerRequirementIndex+i))-*(DS+IndexMaxP+1)+*(DS+i+1))<=DSt)){	//DSt*200 for 1 second
                               dBM[*(SortedPowerRequirementIndex+i)*Ntp][IndexMaxP]=0;
                               *(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i))=(*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i)))-*(P+*(SortedPowerRequirementIndex+i)*Ntp+IndexMaxP)+(*(SortedPowerRequirement+i));
                               *(CumulativeDS+*(SortedPowerRequirementIndex+i))=*(CumulativeDS+*(SortedPowerRequirementIndex+i))-*(DS+IndexMaxP+1)+*(DS+i+1);
                               dBM[*(SortedPowerRequirementIndex+i)*Ntp][i]=1;
                               }
                               else{
                              j=1;
                              while(j<Tbs){
                                            if((*(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i+Ntp*j))+*(SortedPowerRequirement+i+Ntp*j))<=MP&&
                                                *(CumulativeDS+*(SortedPowerRequirementIndex+i+Ntp*j))+*(DS+i+1)<=DSt){	//DSt*200 for 1 second
                                                //printf("Before assignment, MS %d needs DS %d, BS %d: %d\n",i,*(DS+i+1),*(SortedPowerRequirementIndex+i+Ntp*j)+1,*(CumulativeDS+*(SortedPowerRequirementIndex+i+Ntp*j)));
                                                *(UnsortedCumulativePower+*(SortedPowerRequirementIndex+i+Ntp*j))+=*(SortedPowerRequirement+i+Ntp*j);
                                                *(CumulativeDS+*(SortedPowerRequirementIndex+i+Ntp*j))+=*(DS+i+1);
                                                dBM[*(SortedPowerRequirementIndex+i+Ntp*j)][j]=1;
                                                CumulativeServingMSs[*(SortedPowerRequirementIndex+i+Ntp*j)]++;
                                                TotalServedMSs++;
                                                //printf("After assignment, MS %d needs DS %d, BS %d: %d\n",i,*(DS+i+1),*(SortedPowerRequirementIndex+i+Ntp*j)+1,*(CumulativeDS+*(SortedPowerRequirementIndex+i+Ntp*j)));
                                                break;
                                                }
                                            else
                                                j++;
                                            }
                              }
                             }
}*/
       

if(TotalServedMSs<ceil(Ntp*(1-BP))){
	printf("error!!Possible Served MSs:%d. ", TotalServedMSs);
	printf("The required number(%g) of total served MSs(%d) can not be achieved!\n", ceil(Ntp*(1-BP)),Ntp);
	exit(1);
}
       for(i=0;i<Tbs;i++)
	if(BStatus[i]==0)
		*(UnsortedCumulativePower+i)+=SBP;
	else
		*(UnsortedCumulativePower+i)+=ABP;
                           
       for(i=0;i<Tbs;i++){
       //printf("Power from BS %d: %lf\n",i+1,*(UnsortedCumulativePower+i));
       TotalPower+=*(UnsortedCumulativePower+i);
       //printf("DS from BS %d: %d\n",i+1,*(CumulativeDS+i));
       TotalDSs+=*(CumulativeDS+i);
       //printf("MSs served by BS %d: %d\n",i+1,*(CumulativeServingMSs+i));
	//printf("Bandwidth from BS %d: %lf\n",i+1,Bandwidth[i]);
       for(j=0;j<Ntp;j++){
//        printf("dBM%d_%d: %d\n",i+1,j+1,dBM[i][j]);
		if(dBM[i][j]==1)
			NOV+=(BR[j+1]/DS[j+1]);
	}
       }
       
       printf("Total Power %lf\n",TotalPower);
       printf("Total DSs: %d\n",TotalDSs);
       printf("Total served MSs: %d\n",TotalServedMSs);
	printf("Total Bandwidth: %lf\n",TotalBandwidth);
       /*for(j=1;j<=Ntp;j++)
       printf("MS %d needs DSs: %d\n",j,*(DS+j));*/
       
       //dBMprinter(dBM,Ntp,Tbs,DS,P);
       for(i=0;i<Tbs;i++){
       printf("Power from BS %d: %lf\n",i+1,*(UnsortedCumulativePower+i));
       printf("DS from BS %d: %d\n",i+1,*(CumulativeDS+i));
       printf("MSs served by BS %d: %d\n",i+1,*(CumulativeServingMSs+i)); 
	printf("Bandwidth from BS %d: %lf\n\n",i+1,Bandwidth[i]);
       }
       printf("Total Power %lf\n",TotalPower);
       printf("Total DSs: %d\n",TotalDSs);
       printf("Total served MSs: %d\n",TotalServedMSs);
	printf("Total Bandwidth: %lf\n",TotalBandwidth);
	printf("Objective Value: %lf\n",TotalBandwidth-TotalDSs-TotalPower);
	printf("New Objective Value: %lf\n",NOV/TotalPower);
       
       end_tick=clock();
       elapsed=(double)(end_tick-start_tick)/CLOCKS_PER_SEC;
       printf("Running Time: %.9f",elapsed);

	coverage=CoverageFinder(dBM[0],ptrD,Tbs,Ntp);

FILE *HeuristicOut;
if ((HeuristicOut=fopen("HeuristicOut.txt", "w")) == NULL)
printf("\n\nerror!Fail to open file!");
else
printf("\n\nOpen HeuristicOut.txt successfully!\n");

fprintf(HeuristicOut,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
for(i=0;i<Tbs;i++)
fprintf(HeuristicOut," \"BS[%d]\" %d %d %g\n",i+1,Xbs[i+1],Ybs[i+1],*(coverage+i));
fprintf(HeuristicOut,"\n\n");
               for(i=0;i<Tbs;i++){
fprintf(HeuristicOut,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i+1,Xbs[i+1],Ybs[i+1],i+1,Xbs[i+1],Ybs[i+1]);
               for(j=0;j<Ntp;j++)
		if(dBM[i][j]==1){
fprintf(HeuristicOut,"BS[%d] %d %d BS[%d] %d %d\n",i+1,Xbs[i+1],Ybs[i+1],i+1,Xbs[i+1],Ybs[i+1]);
fprintf(HeuristicOut,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW\n",i+1,Xbs[i+1],Ybs[i+1],j+1,Xtp[j+1],Ytp[j+1],DS[j+1],*(P+(i+1)*(Ntp+1)+j+1));
}
fprintf(HeuristicOut,"\n\n");
}
fprintf(HeuristicOut," plot \"HeuristicOut.txt\" index 0:0 using 2:3:1 notitle with labels, \"HeuristicOut.txt\" index 0:0 using 2:3:4 title \"cell size\" with circles, \"coordinates.txt\" index 0:0 using 6:7 title \"users\" with points");
for(i=1;i<=Tbs;i++)
fprintf(HeuristicOut,", \"HeuristicOut.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i,i,i);
fclose(HeuristicOut);
                                  


       return 0;
       }
