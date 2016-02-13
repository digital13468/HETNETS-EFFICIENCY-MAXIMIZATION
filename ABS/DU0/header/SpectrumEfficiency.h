#include<stdio.h> 
#include<time.h>
#include<stdlib.h>
#include"S2Dquicksort4.h"
#include"Squicksort2.h"
#include"3Dquicksort.h"
#include"SpecialQSort.h"
#include"INT_SpecialQSort.h"
#include"quicksort2.h"
#include"ComputingHeader.h"
/*#define SBP 8.00	//Power for BS in sleep mode
#define ABP 68.73	//Basic power for BS in active mode
#define mSBP 7.3	//Sleep mode power level for micro BSs
#define mABP 12	//basic active mode power level for micro BSs*/

/*double * SCoverageFinder(int *dBM, double *D, int Tbs, int Ntp, int Modulation){
	int i,j,k;
	double *coverage=malloc(Tbs*sizeof(double));
	SArrayInitialization(coverage,Tbs);
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(*(dBM+i*Ntp*Modulation+j*Modulation+k)==1)
					if(*(D+i*Ntp+j)/20>coverage[i])
						coverage[i]=*(D+i*Ntp+j)/20;
					//printf("%g\n",coverage[i]);
	return coverage;
}
Stemp_dBM_Initiator(int *ptr, int Ntp){	//initiate the array of temp_dBM for the new heuristic algorithm
	int i,j;
	for(i=0;i<2;i++)
		for(j=0;j<Ntp;j++)
			*(ptr+i*Ntp+j)=-1;
}
SArrayInitialization(double *ptr, int Tbs){
       	int i;         
       	for(i=0;i<Tbs;i++)
       		*(ptr+i)=0.0;
}
SArrayInitialization1(int *ptr, int Tbs){
       	int i;         
       	for(i=0;i<Tbs;i++)
       		*(ptr+i)=0;
}     
SArrayInitialization2(int *ptr, int Tbs, int Ntp, int Modulation){
       	int i,j,k;         
       	for(i=0;i<Tbs;i++)
       		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
       				*(ptr+i*Ntp*Modulation+j*Modulation+k)=0;
}*/           
Sheuristic(int Tbs, double *P, int Ntp, double MP, int DSt, double BP, int *DS, int *Xbs, int *Ybs, int *Xtp, int *Ytp, double *ptrD, double *BR, int Modulation, double ABP, double SBP, double mABP, double mSBP, int mNbs, double mP, int *AdjacentBSs, int *minM, int *Mm, int nUser){
                
       //double *SortedPowerRequirement;  
       //int *SortedPowerRequirementIndex;
	/*double UnsortedObjective[Tbs][Ntp][Modulation];
	SArrayInitialization3(UnsortedObjective[0],Tbs,Ntp,Modulation);
	double SortedObjective[Tbs][Ntp][Modulation];
	int *SortedObjectiveIndex;*/

	double UnsortedCumulativeObjective[Tbs];
	SArrayInitialization(UnsortedCumulativeObjective,Tbs);
// 	double SortedCumulativeObjective[Tbs];
// 	int *SortedCumulativeObjectiveIndex;
       
       	double UnsortedCumulativePower[Tbs];
       	SArrayInitialization(UnsortedCumulativePower,Tbs);
       //double SortedCumulativePower[Tbs];
       	double TotalPower=0.0;
       //int *SortedCumulativePowerIndex;
       
       	int CumulativeDS[Tbs];
       	SArrayInitialization1(CumulativeDS,Tbs);
       	double BSSE=0;

	double CumulativeBR[Tbs];
	SArrayInitialization(CumulativeBR,Tbs);
	double BSP=0.0;
	double NOV=0.0;

	int BStatus[Tbs];
	SArrayInitialization1(BStatus,Tbs);
                            
       	int dBM[Tbs][Ntp][Modulation];
       	int* ptr_dBM=&dBM[0][0][0];
       	SArrayInitialization2(ptr_dBM,Tbs,Ntp,Modulation);
       	int CumulativeServingMSs[Tbs];
       	SArrayInitialization1(CumulativeServingMSs,Tbs);
       	int TotalServedMSs=0;
       
       	int i,j,k,l;
       //double Max=-10000000000.0;
       //int IndexMax=-1;

       	double *coverage;
           
       	clock_t start_tick, end_tick;
       	double elapsed;
       
       	start_tick=clock();

	/*for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++){
			UnsortedObjective[i][j]=(*(BR+(j+1)))/(*(DS+(j+1)))/(*(P+(i+1)*(Ntp+1)+j+1));
			SortedObjective[i][j]=(*(BR+(j+1)))/(*(DS+(j+1)))/(*(P+(i+1)*(Ntp+1)+j+1));
			//printf("BS%d MS%d: %lf \n",i+1,j+1,UnsortedObjective[i][j]);
 			//printf("BR:%lf DS:%d Power:%lf\n",*(BR+(j+1)),*(DS+(j+1)),*(P+(i+1)*(Ntp+1)+j+1));
		}
	//SortedObjective=Unsorted_to_Sorted(UnsortedObjective,Tbs,Ntp);
	SortedObjectiveIndex=Sq_sort_2D(*SortedObjective,0,Tbs-1,Ntp,Tbs);*/
       //SortedPowerRequirement=SortedPowerRequirementAssignment(P,Tbs,Ntp);
       		/*	
                         for(i=0;i<Tbs;i++)
                         for(j=0;j<Ntp;j++)
                         printf("BS%d MS%d: %lf \n",i+1,j+1,SortedObjective[i][j]);
                         */
       //SortedPowerRequirementIndex=q_sort_2D(SortedPowerRequirement,0,Tbs-1,Ntp,Tbs);
	double *SE=malloc(Tbs*Ntp*Modulation*sizeof(double));
	//double SSE[Tbs][Ntp][Modulation];
	//double *SSE=malloc(Tbs*Ntp*Modulation*sizeof(double));
	//double *ptrSSE=&SSE[0][0][0];
	//double *ptrSSE=SSE;
	double *SIndex;
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				
				//SSE[i][j][k]=0;
				*(SE+i*Ntp*Modulation+j*Modulation+k)=0.0;
				//*(SSE+i*Ntp*Modulation+j*Modulation+k)=0;
				//printf("BS%d MS%d: %lf \n",i+1,j+1,UnsortedObjective[i][j]);
				//printf("BR:%lf DS:%d Power:%lf\n",*(BR+(j+1)),*(DS+(j+1)),*(P+(i+1)*(Ntp+1)+j+1));
			}
       for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				*(SE+i*Ntp*Modulation+j*Modulation+k)=((*(BR+(j)))/(*(DS+(j*Modulation)+k)))/(*(P+(i)*(Ntp*Modulation)+j*Modulation+k));
				//*(SSE+i*Ntp*Modulation+j*Modulation+k)=(*(BR+(j)))/(*(DS+(j*Modulation)+k))/(*(P+(i)*(Ntp*Modulation)+j*Modulation+k));
				//SSE[i][j][k]=(*(BR+(j)))/(*(DS+(j*Modulation)+k))/(*(P+(i)*(Ntp*Modulation)+j*Modulation+k));
				//printf(" %lf \n",SE[i][j][k]);
				//printf("BR:%lf DS:%d Power:%lf\n",*(BR+(j+1)),*(DS+(j+1)),*(P+(i+1)*(Ntp+1)+j+1));
			}
      	SIndex=q_sort_3D(SE,Tbs,Ntp,Modulation);
      	int m,n;
      	FILE *poolist;
	if ((poolist=fopen("PooList.txt", "w")) == NULL)
		printf("\n\nerror!Fail to open file!");
	else
		printf("\n\nOpen PooList.txt successfully!\n");
	fprintf(poolist,"#BS%dMS%dBP%gNU%lf\n",Tbs,Ntp,BP, (float)nUser/(float)Ntp);
      	for(n=0;n<Ntp;n++){
		for(m=0;m<Tbs*Modulation;m++){
			fprintf(poolist,"%d MS%d and BS%d use modulation %d transmitting power %lf with se %lf. \n",n*Tbs*Modulation+m,n,(int)floor(*(SIndex+m*Ntp+n)),(int)(*(SIndex+m*Ntp+n)*10-floor(*(SIndex+m*Ntp+n))*10),*(P+(int)floor(*(SIndex+m*Ntp+n))*Ntp*Modulation+n*Modulation+(int)(*(SIndex+m*Ntp+n)*10-floor(*(SIndex+m*Ntp+n))*10)),*(SE+(int)floor(*(SIndex+m*Ntp+n))*Ntp*Modulation+n*Modulation+(int)(*(SIndex+m*Ntp+n)*10-floor(*(SIndex+m*Ntp+n))*10)));
			//printf("%d MS%d %lf transmitting power %lf. \n",n*Tbs*Modulation+m,n,(*(IndexBSMod+n*(Tbs*Modulation)+m)),*(BSMod+n*(Tbs*Modulation)+m));
		}
		fprintf(poolist,"\n");
	}
       	fclose(poolist);
       //spectrum algorithm goes from here
       	int *temp_dBM=malloc(2*Ntp*sizeof(int));
       	Stemp_dBM_Initiator(temp_dBM,Ntp);
// 	double Min=1000.0;
// 	int IndexMin=-1;	
       	for(i=nUser;i<Ntp;i++){	//each MS is assigned to the 1st BS choice
		*(temp_dBM+i)=(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i));
		*(temp_dBM+i+Ntp)=(int)(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)*10-(floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))*10));
		CumulativeServingMSs[(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))]++;
		TotalServedMSs++;
		*(UnsortedCumulativePower+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)))+=*(P+((int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))*Ntp*Modulation)+(i*Modulation)+(int)(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)*10-(floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))*10)));
		*(CumulativeDS+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)))+=*(DS+i*Modulation+(int)(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)*10-(floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))*10)));
		*(CumulativeBR+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)))+=*(BR+i);
		//printf("%lf %lf\n",*(UnsortedCumulativeObjective+*(SortedObjectiveIndex+i+Ntp*(Tbs-1))),SortedObjective[Tbs-1][i]);
		*(UnsortedCumulativeObjective+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)))+=*(SE+((int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))*Ntp*Modulation)+(i*Modulation)+(int)(*(SIndex+Tbs*Ntp*Modulation-Ntp+i)*10-(floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))*10)));
		//printf("%lf\n",*(UnsortedCumulativeObjective+*(SortedObjectiveIndex+i+Ntp*(Tbs-1))));
		BStatus[(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+i))]=1;
	}
	int temp=0;
	for(i=0;i<Tbs;i++){
		printf("%d ", BStatus[i]);
	 //printf("Unsorted Cumulative Objective of BS %d: %lf \n",i+1,*(UnsortedCumulativeObjective+i));
	}
	for(i=0;i<Tbs;i++){
		if(i<Tbs-mNbs&&*(UnsortedCumulativePower+i)>MP){
			printf("Macro BS%d consumes too much energy. Program Stops!",i);
			exit(1);
		}
		if(i>=Tbs-mNbs&&*(UnsortedCumulativePower+i)>mP){
			printf("Micro BS%d consumes too much energy. Program Stops!",i);
			exit(1);
		}
		if(i<Tbs-mNbs){
			temp=*(CumulativeDS+i);
			for(j=Tbs-mNbs;j<Tbs;j++)
				if(*(minM+i*Tbs+j)==1){
					temp+=*(CumulativeDS+j);
						if(temp>DSt){
							printf("Macro BS%d shares too many subcarriers. Program Stops!",i);
							exit(1);
						}
				}
		}
		
	}
	for(j=nUser;j<Ntp;j++)
		if(*(AdjacentBSs+(*(temp_dBM+j))*Ntp+j)==0){
			printf("User%d is out of range of BS%d",j, (*(temp_dBM+j)));
			exit(1);
		}
	int *ServingUserNumber=malloc(Tbs*sizeof(int));
	for(i=0;i<Tbs;i++){
		printf("Total served user %d by BS%d\n",CumulativeServingMSs[i],i);
		*(ServingUserNumber+i)=CumulativeServingMSs[i];
	}
	int *SIndexofUserNumber;
	SIndexofUserNumber=INT_q_sort(ServingUserNumber,0,Tbs-1,Tbs);
// 	for(i=0;i<Tbs;i++)
// 		if(*(SIndexofUserNumber+i)<Tbs-mNbs)
// 			if(TotalServedMSs-*(CumulativeServingMSs+*(SIndexofUserNumber+i))>=ceil(Ntp*(1-BP))){
// 				printf("%d users will be blocked by Macro cell%d.\n",CumulativeServingMSs[*(SIndexofUserNumber+i)],*(SIndexofUserNumber+i));
// 				for(j=0;j<Ntp;j++)
// 					if(*(temp_dBM+j)==*(SIndexofUserNumber+i)){
// 						*(temp_dBM+j)=-1;
// 						CumulativeServingMSs[*(SIndexofUserNumber+i)]--;
// 						TotalServedMSs--;
// 						*(UnsortedCumulativePower+*(SIndexofUserNumber+i))-=*(P+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(CumulativeDS+*(SIndexofUserNumber+i))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
// 						*(CumulativeBR+*(SIndexofUserNumber+i))-=*(BR+j);
// 						*(UnsortedCumulativeObjective+*(SIndexofUserNumber+i))-=*(SE+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(temp_dBM+j+Ntp)=-1;
// 						BStatus[*(SIndexofUserNumber+i)]=0;
// 					}
// 			}
// 	for(i=0;i<Tbs;i++)
// 		if(*(SIndexofUserNumber+i)>=Tbs-mNbs)
// 			if(TotalServedMSs-*(CumulativeServingMSs+*(SIndexofUserNumber+i))>=ceil(Ntp*(1-BP))){
// 				printf("%d users will be blocked by micro cell%d.\n",CumulativeServingMSs[*(SIndexofUserNumber+i)],*(SIndexofUserNumber+i));
// 				for(j=0;j<Ntp;j++)
// 					if(*(temp_dBM+j)==*(SIndexofUserNumber+i)){
// 						*(temp_dBM+j)=-1;
// 						CumulativeServingMSs[*(SIndexofUserNumber+i)]--;
// 						TotalServedMSs--;
// 						*(UnsortedCumulativePower+*(SIndexofUserNumber+i))-=*(P+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(CumulativeDS+*(SIndexofUserNumber+i))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
// 						*(CumulativeBR+*(SIndexofUserNumber+i))-=*(BR+j);
// 						*(UnsortedCumulativeObjective+*(SIndexofUserNumber+i))-=*(SE+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(temp_dBM+j+Ntp)=-1;
// 						BStatus[*(SIndexofUserNumber+i)]=0;
// 					}
// 			}
	int *Mtom=malloc((Tbs-mNbs)*sizeof(int));
	int *MigrableMicro=malloc(Ntp*sizeof(int));
	int TurnOn[Tbs];
	int TotalOn;
	temp=0;
	for(i=Tbs-mNbs;i<Tbs;i++)
		TurnOn[i]=0;
	for(i=0;i<Ntp;i++)
		*(MigrableMicro+i)=-1;
	for(i=0;i<Tbs-mNbs;i++){
		TotalOn=0;
		for(j=Tbs-mNbs;j<Tbs;j++)
			if(*(minM+i*Tbs+j)==1&&BStatus[j]==1)
				TotalOn++;
		*(Mtom+i)=0;
		TurnOn[i]=0;
		if(*(Mm+i)>0){
			for(j=0;j<Ntp;j++)	
				if(*(temp_dBM+j)==i)
					for(k=Tbs-mNbs;k<Tbs;k++)
						if(*(minM+i*Tbs+k)==1&&*(AdjacentBSs+k*Ntp+j)==1){
							temp++;
							*(Mtom+i)=*(Mtom+i)+1;
							if(*(MigrableMicro+j)==-1){
								*(MigrableMicro+j)=k;
								if(TurnOn[k]==0)
									TotalOn++;
								TurnOn[k]++;
							}
							else{
								if((*(P+*(MigrableMicro+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp))>*(P+k*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp)))){
									TurnOn[*(MigrableMicro+j)]--;
									if(TurnOn[*(MigrableMicro+j)]==0)
										TotalOn--;
									*(MigrableMicro+j)=k;
									if(TurnOn[k]==0)
										TotalOn++;
									TurnOn[k]++;
								}	
							}
							
						}
// 			if(TotalServedMSs-(CumulativeServingMSs[i]-temp)<ceil(Ntp*(1-BP))||TotalOn*mABP>ABP){
// 				temp-=*(Mtom+i);
// 				*(Mtom+i)=0;
// 			}
			if(*(Mtom+i)!=CumulativeServingMSs[i]||TotalOn*mABP>ABP){
				temp-=*(Mtom+i);
				*(Mtom+i)=0;
			}
		}
	}
	for(i=0;i<Tbs-mNbs;i++)
		printf("BS%d has %d user to be moved.\n",i,*(Mtom+i));
	printf("The total movable user number is %d.\n",temp);
	if(temp>0){
		int *SortedMtom;
		int *ModulationSelection=malloc(Ntp*sizeof(int));
		int *TotalSubcarrier=malloc(Tbs*sizeof(int));
		double *TotalEnergy=malloc(Tbs*sizeof(double));
		int *temp_BStatus=malloc(Tbs*sizeof(int));
		int *mBelongsM=malloc(Tbs*sizeof(int));
		int *Migrato=malloc(Ntp*sizeof(int));
		int *ServedUser=malloc(Tbs*sizeof(int));
		int TrueMigration=0;
		for(j=0;j<Ntp;j++){
			*(ModulationSelection+j)=*(temp_dBM+j+Ntp);
			*(Migrato+j)=*(temp_dBM+j);
		}
		for(j=0;j<Tbs;j++){
			*(ServedUser+j)=CumulativeServingMSs[j];
			*(TotalEnergy+j)=*(UnsortedCumulativePower+j);
			*(temp_BStatus+j)=BStatus[j];
			*(TotalSubcarrier+j)=*(CumulativeDS+j);
			if(j>=Tbs-mNbs){
				for(k=0;k<Tbs-mNbs;k++)
					if(*(minM+k*Tbs+j)==1){
						*(TotalSubcarrier+k)+=*(TotalSubcarrier+j);
						*(mBelongsM+j)=k;
					}
			}
			else
				*(mBelongsM+j)=-1;
		}
		int MigrableUser;
		SortedMtom=INT_q_sort(Mtom,0,Tbs-mNbs-1,Tbs-mNbs);
		for(i=0;i<Tbs-mNbs;i++){
			MigrableUser=0;
			if(*(Mtom+i)>0){
				for(j=0;j<Ntp;j++)
					if(*(minM+*(SortedMtom+i)*Tbs+*(MigrableMicro+j))==1&&*(temp_dBM+j)==*(SortedMtom+i)){
						if(*(TotalEnergy+*(MigrableMicro+j))+*(P+*(MigrableMicro+j)*Ntp*Modulation+j*Modulation+*(ModulationSelection+j))<=mP){
							MigrableUser++;
							*(Migrato+j)=*(MigrableMicro+j);
							*(TotalEnergy+*(MigrableMicro+j))+=*(P+*(MigrableMicro+j)*Ntp*Modulation+j*Modulation+*(ModulationSelection+j));
							*(ServedUser+*(MigrableMicro+j))=*(ServedUser+*(MigrableMicro+j))+1;
							*(TotalSubcarrier+*(MigrableMicro+j))+=*(DS+j*Modulation+*(ModulationSelection+j));
							*(TotalEnergy+*(SortedMtom+i))-=*(P+*(SortedMtom+i)*Ntp*Modulation+j*Modulation+*(ModulationSelection+j));
							*(ServedUser+*(SortedMtom+i))=*(ServedUser+*(SortedMtom+i))-1;
							if(*(temp_BStatus+*(MigrableMicro+j))==0)
								*(temp_BStatus+*(MigrableMicro+j))=1;
						}	
						else{
							*(Migrato+j)=*(MigrableMicro+j);
							*(ServedUser+*(MigrableMicro+j))=*(ServedUser+*(MigrableMicro+j))+1;
							*(TotalEnergy+*(MigrableMicro+j))+=*(P+*(MigrableMicro+j)*Ntp*Modulation+j*Modulation+*(ModulationSelection+j));
							*(TotalSubcarrier+*(MigrableMicro+j))+=*(DS+j*Modulation+*(ModulationSelection+j));
							*(TotalEnergy+*(SortedMtom+i))-=*(P+*(SortedMtom+i)*Ntp*Modulation+j*Modulation+*(ModulationSelection+j));
							*(ServedUser+*(SortedMtom+i))=*(ServedUser+*(SortedMtom+i))-1;
							double *EnergyCon=malloc(*(ServedUser+*(MigrableMicro+j))*sizeof(double));
							int *EnergyConIndex=malloc(*(ServedUser+*(MigrableMicro+j))*sizeof(int));
							l=0;
							for(k=0;k<Ntp;k++)
								if(*(Migrato+k)==*(MigrableMicro+j)){
									*(EnergyCon+l)=*(P+*(MigrableMicro+j)*Ntp*Modulation+k*Modulation+*(ModulationSelection+k));
									*(EnergyConIndex+l)=k;
									l++;
								}
							DoubleIndex_q_sort(EnergyCon,0,*(ServedUser+*(MigrableMicro+j))-1,*(ServedUser+*(MigrableMicro+j)),EnergyConIndex);
							l=*(ServedUser+*(MigrableMicro+j))-1;
							do{
								if(*(ModulationSelection+*(EnergyConIndex+l))==0)
									l--;
								else{
									*(ModulationSelection+*(EnergyConIndex+l))=*(ModulationSelection+*(EnergyConIndex+l))-1;
									*(TotalEnergy+*(MigrableMicro+j))=*(TotalEnergy+*(MigrableMicro+j))+*(P+*(MigrableMicro+j)*Ntp*Modulation+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)))-*(P+*(MigrableMicro+j)*Ntp*Modulation+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l))+1);
									*(TotalSubcarrier+*(SortedMtom+i))=*(TotalSubcarrier+*(SortedMtom+i))+*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)))-*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l))+1);
									*(TotalSubcarrier+*(MigrableMicro+j))=*(TotalSubcarrier+*(MigrableMicro+j))+*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)))-*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l))+1);
								}
							}while(*(TotalEnergy+*(MigrableMicro+j))>mP&&l>-1&&*(TotalSubcarrier+*(SortedMtom+i))<=DSt);
							if(*(TotalEnergy+*(MigrableMicro+j))<=mP&&l>-1&&*(TotalSubcarrier+*(SortedMtom+i))<=DSt)
								MigrableUser++;
							else{
								*(Migrato+j)=*(temp_dBM+j);
								*(ServedUser+*(MigrableMicro+j))=*(ServedUser+*(MigrableMicro+j))-1;
								*(TotalEnergy+*(MigrableMicro+j))-=*(P+*(MigrableMicro+j)*Ntp*Modulation+j*Modulation+*(ModulationSelection+j));
								*(TotalSubcarrier+*(MigrableMicro+j))-=*(DS+j*Modulation+*(ModulationSelection+j));
								*(TotalEnergy+*(SortedMtom+i))+=*(P+*(SortedMtom+i)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
								*(TotalSubcarrier+*(SortedMtom+i))=*(TotalSubcarrier+*(SortedMtom+i))+*(DS+j*Modulation+*(temp_dBM+j+Ntp))-*(DS+j*Modulation+*(ModulationSelection+j));
								*(ServedUser+*(SortedMtom+i))=*(ServedUser+*(SortedMtom+i))+1;
								*(ModulationSelection+j)=*(temp_dBM+j+Ntp);
								if(*(TotalSubcarrier+*(SortedMtom+i))>DSt&&*(EnergyConIndex+l)!=j){
									*(ModulationSelection+*(EnergyConIndex+l))=*(ModulationSelection+*(EnergyConIndex+l))+1;
									*(TotalEnergy+*(MigrableMicro+j))=*(TotalEnergy+*(MigrableMicro+j))-*(P+*(MigrableMicro+j)*Ntp*Modulation+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)-1))+*(P+*(MigrableMicro+j)*Ntp*Modulation+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)));
									*(TotalSubcarrier+*(MigrableMicro+j))=*(TotalSubcarrier+*(MigrableMicro+j))-*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)-1))+*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)));
									*(TotalSubcarrier+*(SortedMtom+i))=*(TotalSubcarrier+*(SortedMtom+i))-*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)-1))+*(DS+*(EnergyConIndex+l)*Modulation+*(ModulationSelection+*(EnergyConIndex+l)));
								}
							}	
							free(EnergyCon);
							free(EnergyConIndex);
						}
					}	
			}
//			if(TotalServedMSs-(CumulativeServingMSs[*(SortedMtom+i)]-MigrableUser)>=ceil(Ntp*(1-BP))){
			if(MigrableUser==CumulativeServingMSs[*(SortedMtom+i)]){
				TrueMigration=1;
				TotalServedMSs-=(CumulativeServingMSs[*(SortedMtom+i)]-MigrableUser);
				for(j=0;j<Ntp;j++){
					if(*(Migrato+j)==*(SortedMtom+i)&&*(Migrato+j)==*(temp_dBM+j)){
						*(temp_dBM+j+Ntp)=-1;
						*(temp_dBM+j)=-1;
					}
					else if(*(minM+*(temp_dBM+j)*Tbs+*(Migrato+j))==1&&*(temp_dBM+j)==*(SortedMtom+i)){
						*(temp_dBM+j+Ntp)=*(ModulationSelection+j);
						*(temp_dBM+j)=*(Migrato+j);
					}
					else if(*(temp_dBM+j)==*(Migrato+j)&&*(minM+*(SortedMtom+i)*Tbs+*(Migrato+j))==1){
						*(temp_dBM+j+Ntp)=*(ModulationSelection+j);
						*(temp_dBM+j)=*(Migrato+j);
					}
				}
				for(j=0;j<Tbs;j++){
					if(*(SortedMtom+i)==j){
						CumulativeServingMSs[j]=0;
						*(UnsortedCumulativePower)=0;
						BStatus[j]=0;
						*(CumulativeDS+j)=0;
						*(TotalSubcarrier+j)=0;
					}
					else if(*(minM+*(SortedMtom+i)*Tbs+j)==1){
						CumulativeServingMSs[j]=*(ServedUser+j);
						*(UnsortedCumulativePower+j)=*(TotalEnergy+j);
						BStatus[j]=*(temp_BStatus+j);
						*(CumulativeDS+j)=*(TotalSubcarrier+j);
						*(TotalSubcarrier+*(SortedMtom+i))+=*(TotalSubcarrier+j);
					}
				}
				if(TotalServedMSs==ceil(Ntp*(1-BP)))
					i=Tbs-mNbs;
			}	
			
		}
		free(Migrato);
		SortedMtom = NULL;
		free(mBelongsM);
		free(ModulationSelection);
		free(temp_BStatus);
		free(ServedUser);
		if(TrueMigration==1){
			int TUC=0;
			int BUC[Tbs];
			SArrayInitialization1(BUC,Tbs);
			double PC[Tbs];
			SArrayInitialization(PC,Tbs);
			int SC[Tbs];
			SArrayInitialization1(SC,Tbs);
			double SEC[Tbs];
			SArrayInitialization(SEC,Tbs);
			int SI[Tbs];
			SArrayInitialization1(SI,Tbs);
			double BC[Tbs];
			SArrayInitialization(BC,Tbs);
			for(i=0;i<Ntp;i++)
				if(*(temp_dBM+i)!=-1){
					TUC++;
					if(*(temp_dBM+i+Ntp)==-1){
						printf("error! User%d has connection but no modulation.\n",i);
						exit(1);
					}
					else if(*(AdjacentBSs+*(temp_dBM+i)*Ntp+i)!=1){
						printf("error! user%d cannot connect with BS%d.\n",i,*(temp_dBM+i));
						exit(1);
					}
					BUC[*(temp_dBM+i)]++;
					PC[*(temp_dBM+i)]=*(P+*(temp_dBM+i)*Ntp*Modulation+i*Modulation+*(temp_dBM+i+Ntp));
					SC[*(temp_dBM+i)]=*(DS+i*Modulation+*(temp_dBM+i+Ntp));
					SEC[*(temp_dBM+i)]=*(SE+*(temp_dBM+i)*Ntp*Modulation+i*Modulation+*(temp_dBM+i+Ntp));
					SI[*(temp_dBM+i)]=1;
					BC[*(temp_dBM+i)]+=*(BR+i);
				}
			if(TUC!=TotalServedMSs){
				printf("error! The numbers of served users in total are different!(TUC%d TotalServedMSs%d)\n",TUC,TotalServedMSs);
				exit(1);
			}
			for(i=0;i<Tbs;i++){
				if(*(TotalSubcarrier+i)>DSt){
					printf("BS%d violates DS constraint\n",i);
					exit(1);
				}
				else if(i<Tbs-mNbs&&*(TotalEnergy+i)>MP){
					printf("BS%d violates power constraint\n",i);
					exit(1);
				}
				else if(i>=Tbs-mNbs&&*(TotalEnergy+i)>mP){
					printf("BS%d violates power constraint\n",i);
					exit(1);
				}
				else if(BUC[i]!=CumulativeServingMSs[i]){
					printf("error! The numbers of served users by BS%d are different!(BUC%d CumulativeServingMSs%d)\n",i,BUC[i],CumulativeServingMSs[i]);
					exit(1);
				}
				else if(PC[i]!=*(UnsortedCumulativePower+i)){
					printf("1. error! The energy consumed by BS%d are different!(PC%g UnsortedCumulativePower%g)\n",i,PC[i],*(UnsortedCumulativePower+i));
					exit(1);
				}
				else if(SC[i]!=*(CumulativeDS+i)){
					printf("error! The subcarrier used by BS%d are different!(SC%d CumulativeDS%d)\n",i,SC[i],*(CumulativeDS+i));
					exit(1);
				}
				
				else if(SI[i]!=BStatus[i]){
					printf("error! The modes of BS%d are different!(SI%d BStatus%d)\n",i,SI[i],BStatus[i]);
					exit(1);
				}
				*(CumulativeBR+i)=BC[i];
				*(UnsortedCumulativeObjective+i)=SEC[i];
			}
		}
		free(TotalSubcarrier);
		free(TotalEnergy);
	}
	free(Mtom);
	free(MigrableMicro);
	int CoverableUsers[Tbs];
	int Remaining=0;
	for(i=0;i<Tbs;i++){
		CoverableUsers[i]=0;
		for(j=0;j<Ntp;j++)
			if(*(AdjacentBSs+i*Ntp+j)==1&&*(temp_dBM+j)!=-1)
				CoverableUsers[i]++;
		if(CoverableUsers[i]>0)
			Remaining++;
	}
	int* SCoverableUsersIndex;
	i=Remaining-1;
	int BSBreath[Tbs];
	SArrayInitialization1(BSBreath,Tbs);
	int MigrateUser[Ntp];
	SArrayInitialization1(MigrateUser,Ntp);
	while(i>-1){
		int* CovNum=malloc((i+1)*sizeof(int));
		int* CovInd=malloc((i+1)*sizeof(int));
		n=0;
		for(m=0;m<i+1;m++){
			if(CoverableUsers[n]>0){
				*(CovNum+m)=CoverableUsers[n];
				*(CovInd+m)=n;
			}
			else
				m--;
			n++;
		}
		SCoverableUsersIndex=Index_q_sort(CovNum,0,i,i+1,CovInd);
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==*(SCoverableUsersIndex+i))
				MigrateUser[j]=1;
			else if(*(AdjacentBSs+*(SCoverableUsersIndex+i)*Ntp+j)==1&&*(temp_dBM+j)!=-1&&*(temp_dBM+j)!=*(SCoverableUsersIndex+i)){
				double *UnsortedUserSE=malloc(Modulation*sizeof(double));
				int* SUserSE;
				for(k=0;k<Modulation;k++)
					*(UnsortedUserSE+k)=(*(BR+(j)))/(*(DS+(j*Modulation)+k))/(*(P+(*(SCoverableUsersIndex+i)*Ntp*Modulation)+j*Modulation+k));
				SUserSE=DOUBLE_q_sort(UnsortedUserSE,0,Modulation-1,Modulation);
					if(*(SCoverableUsersIndex+i)<Tbs-mNbs){
						temp=0;
						for(l=Tbs-mNbs;l<Tbs;l++)
							if(*(minM+*(SCoverableUsersIndex+i)*Tbs+l)==1)
								temp+=*(CumulativeDS+l);
						temp+=*(CumulativeDS+*(SCoverableUsersIndex+i));
						if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+(*(SCoverableUsersIndex+i)*Ntp*Modulation)+j*Modulation+*(SUserSE+Modulation-1))<=MP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))<=DSt){
							CumulativeServingMSs[*(temp_dBM+j)]--;
							TotalServedMSs--;
							*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
							*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
							*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
							*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
							if(CumulativeServingMSs[*(temp_dBM+j)]==0)
								BStatus[*(temp_dBM+j)]=0;
							*(temp_dBM+j)=*(SCoverableUsersIndex+i);
							*(temp_dBM+j+Ntp)=*(SUserSE+Modulation-1);
							CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
							TotalServedMSs++;
							*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
							*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
							*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
							*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
							BStatus[*(SCoverableUsersIndex+i)]=1;
							MigrateUser[j]=1;
						}
						else if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+(*(SCoverableUsersIndex+i)*Ntp*Modulation)+j*Modulation+*(SUserSE+Modulation-1))>MP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))<=DSt){
							int TuneModulation=*(SUserSE+Modulation-1);
							while(TuneModulation-1>=0){
								if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+(*(SCoverableUsersIndex+i)*Ntp*Modulation)+j*Modulation+TuneModulation-1)<=MP&&temp+*(DS+j*Modulation+TuneModulation-1)<=DSt){
									CumulativeServingMSs[*(temp_dBM+j)]--;
									TotalServedMSs--;
									*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
									*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									if(CumulativeServingMSs[*(temp_dBM+j)]==0)
										BStatus[*(temp_dBM+j)]=0;
									*(temp_dBM+j)=*(SCoverableUsersIndex+i);
									*(temp_dBM+j+Ntp)=TuneModulation-1;
									CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
									TotalServedMSs++;
									*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1);
									*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+TuneModulation-1);
									*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
									*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1);
									BStatus[*(SCoverableUsersIndex+i)]=1;
									MigrateUser[j]=1;
									TuneModulation=-1;
								}
								else if(TuneModulation==1&&*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1)>MP&&temp+*(DS+j*Modulation+TuneModulation-1)<DSt){
									double *ExistingUserP=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(double));
									int *ExistingUserPIndex=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(int));
									int x,y=0;
									for(x=0;x<Ntp;x++)
										if(*(temp_dBM+x)==*(SCoverableUsersIndex+i)){
											*(ExistingUserP+y)=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+x*Modulation+*(temp_dBM+x+Ntp));
											*(ExistingUserPIndex+y)=x;
											y++;
										}
									*(ExistingUserP+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
									*(ExistingUserPIndex+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=j;
									int *SExistingUserP;
									SExistingUserP=DoubleIndex_q_sort(ExistingUserP,0,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1-1,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1,ExistingUserPIndex);
									int Tuning=0;
									int z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
									double TPower=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
									double TDS=temp+*(DS+j*Modulation+*(SUserSE+Modulation-1));
									while(Tuning==0&&z>=0){
											int TMod=*(temp_dBM+*(SExistingUserP+z)+Ntp);
											while(TMod>0&&Tuning==0){
												TPower=TPower-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+TMod)+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+TMod-1);
												TDS=TDS-*(DS+*(SExistingUserP+z)*Modulation+TMod)+*(DS+*(SExistingUserP+z)*Modulation+TMod-1);
												if(TPower<=MP&&TDS<=DSt)
													Tuning=1;
												else
													TMod--;	
											}
										z--;
									}
									if(Tuning==1){
										CumulativeServingMSs[*(temp_dBM+j)]--;
										TotalServedMSs--;
										*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
										*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										if(CumulativeServingMSs[*(temp_dBM+j)]==0)
											BStatus[*(temp_dBM+j)]=0;
										*(temp_dBM+j)=*(SCoverableUsersIndex+i);
										*(temp_dBM+j+Ntp)=*(SUserSE+Modulation-1);
										CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
										TotalServedMSs++;
										*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
										*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										BStatus[*(SCoverableUsersIndex+i)]=1;
										MigrateUser[j]=1;
										temp+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
										while(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))>MP||temp>DSt){
											if(*(temp_dBM+*(SExistingUserP+z)+Ntp)>0){
												*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												*(CumulativeDS+*(SCoverableUsersIndex+i))=*(CumulativeDS+*(SCoverableUsersIndex+i))-*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))=*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))-*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												temp=temp-*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												*(temp_dBM+*(SExistingUserP+z)+Ntp)=*(temp_dBM+*(SExistingUserP+z)+Ntp)-1;
											}
											else
												z--;
										}
									}
									free(ExistingUserP);
									free(ExistingUserPIndex);
									TuneModulation=-1;
								}
								else
									TuneModulation-=1;
							}
						}	
						else if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1))<=MP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))>DSt){
							int TuneModulation=*(SUserSE+Modulation-1);
							while(TuneModulation+1<Modulation){
								if(temp+*(DS+j*Modulation+TuneModulation+1)<=DSt&&*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation+1)<=MP){
									CumulativeServingMSs[*(temp_dBM+j)]--;
									TotalServedMSs--;
									*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
									*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									if(CumulativeServingMSs[*(temp_dBM+j)]==0)
										BStatus[*(temp_dBM+j)]=0;
									*(temp_dBM+j)=*(SCoverableUsersIndex+i);
									*(temp_dBM+j+Ntp)=TuneModulation+1;
									CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
									TotalServedMSs++;
									*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation+1);
									*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+TuneModulation+1);
									*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
									*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation+1);
									BStatus[*(SCoverableUsersIndex+i)]=1;
									MigrateUser[j]=1;	
									TuneModulation=1000;
								}
								else if(TuneModulation==Modulation-2&&*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1)<MP&&temp+*(DS+j*Modulation+TuneModulation-1)>DSt){
									int *ExistingUserD=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(int));
									int *ExistingUserDIndex=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(int));
									int x,y=0;
									for(x=0;x<Ntp;x++)
										if(*(temp_dBM+x)==*(SCoverableUsersIndex+i)){
											*(ExistingUserD+y)=*(DS+x*Modulation+*(temp_dBM+x+Ntp));
											*(ExistingUserDIndex+y)=x;
											y++;
										}
									*(ExistingUserD+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=*(DS+j*Modulation+*(SUserSE+Modulation-1));
									*(ExistingUserDIndex+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=j;
									int *SExistingUserD;
									SExistingUserD=Index_q_sort(ExistingUserD,0,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1-1,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1,ExistingUserDIndex);
									int Tuning=0;
									int z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
									double TPower=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
									double TDS=temp+*(DS+j*Modulation+*(SUserSE+Modulation-1));
									while(Tuning==0&&z>=0){
										int TMod=*(temp_dBM+*(SExistingUserD+z)+Ntp);
										while(TMod<Modulation-1&&Tuning==0){
											TPower=TPower-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+TMod)+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+TMod+1);
											TDS=TDS-*(DS+*(SExistingUserD+z)*Modulation+TMod)+*(DS+*(SExistingUserD+z)*Modulation+TMod+1);
											if(TPower<=MP&&TDS<=DSt)
												Tuning=1;
											else
												TMod++;	
										}
										z--;
									}
									if(Tuning==1){
										CumulativeServingMSs[*(temp_dBM+j)]--;
										TotalServedMSs--;
										*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
										*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										if(CumulativeServingMSs[*(temp_dBM+j)]==0)
											BStatus[*(temp_dBM+j)]=0;
										*(temp_dBM+j)=*(SCoverableUsersIndex+i);
										*(temp_dBM+j+Ntp)=*(SUserSE+Modulation-1);
										CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
										TotalServedMSs++;
										*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
										*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										BStatus[*(SCoverableUsersIndex+i)]=1;
										MigrateUser[j]=1;
										temp+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
										while(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))>MP||temp>DSt){
											if(*(temp_dBM+*(SExistingUserD+z)+Ntp)<Modulation-1){
												*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												*(CumulativeDS+*(SCoverableUsersIndex+i))=*(CumulativeDS+*(SCoverableUsersIndex+i))-*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))=*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))-*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												temp=temp-*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												*(temp_dBM+*(SExistingUserD+z)+Ntp)=*(temp_dBM+*(SExistingUserD+z)+Ntp)+1;
											}
											else
												z--;
										}
									}
									free(ExistingUserD);
									free(ExistingUserDIndex);
									TuneModulation=1000;
								}
								else
									TuneModulation+=1;
							}
						}
						//else if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1))>MP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))>DSt)
					}
					else if(*(SCoverableUsersIndex+i)>=Tbs-mNbs){
						temp=0;
						int CoverMBS=-1;
						for(l=0;l<Tbs-mNbs;l++)
							if(*(minM+l*Tbs+*(SCoverableUsersIndex+i))==1)
								CoverMBS=l;
						for(l=Tbs-mNbs;l<Tbs;l++)
							if(*(minM+CoverMBS*Tbs+l)==1)
								temp+=*(CumulativeDS+l);
						temp+=*(CumulativeDS+CoverMBS);
						if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+(*(SCoverableUsersIndex+i)*Ntp*Modulation)+j*Modulation+*(SUserSE+Modulation-1))<=mP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))<=DSt){
							CumulativeServingMSs[*(temp_dBM+j)]--;
							TotalServedMSs--;
							*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
							*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
							*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
							*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
							if(CumulativeServingMSs[*(temp_dBM+j)]==0)
								BStatus[*(temp_dBM+j)]=0;
							*(temp_dBM+j)=*(SCoverableUsersIndex+i);
							*(temp_dBM+j+Ntp)=*(SUserSE+Modulation-1);
							CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
							TotalServedMSs++;
							*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
							*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
							*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
							*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
							BStatus[*(SCoverableUsersIndex+i)]=1;
							MigrateUser[j]=1;
						}
						else if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+(*(SCoverableUsersIndex+i)*Ntp*Modulation)+j*Modulation+*(SUserSE+Modulation-1))>mP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))<=DSt){
							int TuneModulation=*(SUserSE+Modulation-1);
							while(TuneModulation-1>=0){
								if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+(*(SCoverableUsersIndex+i)*Ntp*Modulation)+j*Modulation+TuneModulation-1)<=mP&&temp+*(DS+j*Modulation+TuneModulation-1)<=DSt){
									CumulativeServingMSs[*(temp_dBM+j)]--;
									TotalServedMSs--;
									*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
									*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									if(CumulativeServingMSs[*(temp_dBM+j)]==0)
										BStatus[*(temp_dBM+j)]=0;
									*(temp_dBM+j)=*(SCoverableUsersIndex+i);
									*(temp_dBM+j+Ntp)=TuneModulation-1;
									CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
									TotalServedMSs++;
									*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1);
									*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+TuneModulation-1);
									*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
									*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1);
									BStatus[*(SCoverableUsersIndex+i)]=1;
									MigrateUser[j]=1;
									TuneModulation=-1;
								}
								else if(TuneModulation==1&&*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1)>mP&&temp+*(DS+j*Modulation+TuneModulation-1)<DSt){
									double *ExistingUserP=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(double));
									int *ExistingUserPIndex=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(int));
									int x,y=0;
									for(x=0;x<Ntp;x++)
										if(*(temp_dBM+x)==*(SCoverableUsersIndex+i)){
											*(ExistingUserP+y)=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+x*Modulation+*(temp_dBM+x+Ntp));
											*(ExistingUserPIndex+y)=x;
											y++;
										}
									*(ExistingUserP+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
									*(ExistingUserPIndex+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=j;
									int *SExistingUserP;
									SExistingUserP=DoubleIndex_q_sort(ExistingUserP,0,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1-1,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1,ExistingUserPIndex);
									int Tuning=0;
									int z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
									double TPower=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
									double TDS=temp+*(DS+j*Modulation+*(SUserSE+Modulation-1));
									while(Tuning==0&&z>=0){
											int TMod=*(temp_dBM+*(SExistingUserP+z)+Ntp);
											while(TMod>0&&Tuning==0){
												TPower=TPower-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+TMod)+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+TMod-1);
												TDS=TDS-*(DS+*(SExistingUserP+z)*Modulation+TMod)+*(DS+*(SExistingUserP+z)*Modulation+TMod-1);
												if(TPower<=mP&&TDS<=DSt)
													Tuning=1;
												else
													TMod--;	
											}
										z--;
									}
									if(Tuning==1){
										CumulativeServingMSs[*(temp_dBM+j)]--;
										TotalServedMSs--;
										*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
										*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										if(CumulativeServingMSs[*(temp_dBM+j)]==0)
											BStatus[*(temp_dBM+j)]=0;
										*(temp_dBM+j)=*(SCoverableUsersIndex+i);
										*(temp_dBM+j+Ntp)=*(SUserSE+Modulation-1);
										CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
										TotalServedMSs++;
										*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
										*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										BStatus[*(SCoverableUsersIndex+i)]=1;
										MigrateUser[j]=1;
										temp+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
										while(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))>mP||temp>DSt){
											if(*(temp_dBM+*(SExistingUserP+z)+Ntp)>0){
												*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												*(CumulativeDS+*(SCoverableUsersIndex+i))=*(CumulativeDS+*(SCoverableUsersIndex+i))-*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))=*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))-*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												temp=temp-*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp))+*(DS+*(SExistingUserP+z)*Modulation+*(temp_dBM+*(SExistingUserP+z)+Ntp)-1);
												*(temp_dBM+*(SExistingUserP+z)+Ntp)=*(temp_dBM+*(SExistingUserP+z)+Ntp)-1;
											}
											else
												z--;
										}
									}
									free(ExistingUserP);
									free(ExistingUserPIndex);
									TuneModulation=-1;
								}
								else
									TuneModulation-=1;
							}
						}	
						else if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1))<=mP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))>DSt){
							int TuneModulation=*(SUserSE+Modulation-1);
							while(TuneModulation+1<Modulation){
								if(temp+*(DS+j*Modulation+TuneModulation+1)<=DSt&&*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation+1)<=mP){
									CumulativeServingMSs[*(temp_dBM+j)]--;
									TotalServedMSs--;
									*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
									*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
									*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
									if(CumulativeServingMSs[*(temp_dBM+j)]==0)
										BStatus[*(temp_dBM+j)]=0;
									*(temp_dBM+j)=*(SCoverableUsersIndex+i);
									*(temp_dBM+j+Ntp)=TuneModulation+1;
									CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
									TotalServedMSs++;
									*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation+1);
									*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+TuneModulation+1);
									*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
									*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation+1);
									BStatus[*(SCoverableUsersIndex+i)]=1;
									MigrateUser[j]=1;	
									TuneModulation=1000;
								}
								else if(TuneModulation==Modulation-2&&*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+TuneModulation-1)<mP&&temp+*(DS+j*Modulation+TuneModulation-1)>DSt){
									int *ExistingUserD=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(int));
									int *ExistingUserDIndex=malloc((CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)*sizeof(int));
									int x,y=0;
									for(x=0;x<Ntp;x++)
										if(*(temp_dBM+x)==*(SCoverableUsersIndex+i)){
											*(ExistingUserD+y)=*(DS+x*Modulation+*(temp_dBM+x+Ntp));
											*(ExistingUserDIndex+y)=x;
											y++;
										}
									*(ExistingUserD+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=*(DS+j*Modulation+*(SUserSE+Modulation-1));
									*(ExistingUserDIndex+(CumulativeServingMSs[*(SCoverableUsersIndex+i)]))=j;
									int *SExistingUserD;
									SExistingUserD=Index_q_sort(ExistingUserD,0,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1-1,CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1,ExistingUserDIndex);
									int Tuning=0;
									int z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
									double TPower=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
									double TDS=temp+*(DS+j*Modulation+*(SUserSE+Modulation-1));
									while(Tuning==0||z<(CumulativeServingMSs[*(SCoverableUsersIndex+i)]+1)){
										int TMod=*(temp_dBM+*(SExistingUserD+z)+Ntp);
										while(TMod<Modulation-1&&Tuning==0){
											TPower=TPower-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+TMod)+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+TMod+1);
											TDS=TDS-*(DS+*(SExistingUserD+z)*Modulation+TMod)+*(DS+*(SExistingUserD+z)*Modulation+TMod+1);
											if(TPower<=mP&&TDS<=DSt)
												Tuning=1;
											else
												TMod++;	
										}
										z--;
									}
									if(Tuning==1){
										CumulativeServingMSs[*(temp_dBM+j)]--;
										TotalServedMSs--;
										*(UnsortedCumulativePower+*(temp_dBM+j))-=*(P+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeDS+*(temp_dBM+j))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
										*(CumulativeBR+*(temp_dBM+j))-=*(BR+j);
										*(UnsortedCumulativeObjective+*(temp_dBM+j))-=*(SE+*(temp_dBM+j)*Ntp*Modulation+j*Modulation+*(temp_dBM+j+Ntp));
										if(CumulativeServingMSs[*(temp_dBM+j)]==0)
											BStatus[*(temp_dBM+j)]=0;
										*(temp_dBM+j)=*(SCoverableUsersIndex+i);
										*(temp_dBM+j+Ntp)=*(SUserSE+Modulation-1);
										CumulativeServingMSs[*(SCoverableUsersIndex+i)]++;
										TotalServedMSs++;
										*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+=*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeDS+*(SCoverableUsersIndex+i))+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										*(CumulativeBR+*(SCoverableUsersIndex+i))+=*(BR+j);
										*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))+=*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1));
										BStatus[*(SCoverableUsersIndex+i)]=1;
										MigrateUser[j]=1;
										temp+=*(DS+j*Modulation+*(SUserSE+Modulation-1));
										z=CumulativeServingMSs[*(SCoverableUsersIndex+i)];
										while(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))>mP||temp>DSt){
											if(*(temp_dBM+*(SExistingUserD+z)+Ntp)<Modulation-1){
												*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))=*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))-*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												*(CumulativeDS+*(SCoverableUsersIndex+i))=*(CumulativeDS+*(SCoverableUsersIndex+i))-*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))=*(UnsortedCumulativeObjective+*(SCoverableUsersIndex+i))-*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(SE+*(SCoverableUsersIndex+i)*Ntp*Modulation+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												temp=temp-*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp))+*(DS+*(SExistingUserD+z)*Modulation+*(temp_dBM+*(SExistingUserD+z)+Ntp)+1);
												*(temp_dBM+*(SExistingUserD+z)+Ntp)=*(temp_dBM+*(SExistingUserD+z)+Ntp)+1;
											}
											else
												z--;
										}
									}
									free(ExistingUserD);
									free(ExistingUserDIndex);
									TuneModulation=1000;
								}
								else
									TuneModulation+=1;
							}
						}
						//else if(*(UnsortedCumulativePower+*(SCoverableUsersIndex+i))+*(P+*(SCoverableUsersIndex+i)*Ntp*Modulation+j*Modulation+*(SUserSE+Modulation-1))>MP&&temp+*(DS+j*Modulation+*(SUserSE+Modulation-1))>DSt)
					}
				free(UnsortedUserSE);
			}
		BSBreath[*(SCoverableUsersIndex+i)]=1;
		Remaining=0;
		for(m=0;m<Tbs;m++){
			CoverableUsers[m]=0;
			for(n=0;n<Ntp;n++)
				if(MigrateUser[n]==0&&BSBreath[m]==0&&*(AdjacentBSs+m*Ntp+n)==1&&*(temp_dBM+n)!=-1)
					CoverableUsers[m]++;
			if(CoverableUsers[m]>0)
				Remaining++;
		}
		if(Remaining>0)
			i=Remaining-1;
		else
			i=-1;
		free(CovNum);
		free(CovInd);
	}
// 	for(i=0;i<Tbs;i++){
// 		printf("Total served user %d by BS%d\n",CumulativeServingMSs[i],i);
// 		*(ServingUserNumber+i)=CumulativeServingMSs[i];
// 	}
// 	SIndexofUserNumber=INT_q_sort(ServingUserNumber,0,Tbs-1,Tbs);
// 	
// 	for(i=0;i<Tbs;i++)
// 		if(*(SIndexofUserNumber+i)<Tbs-mNbs)
// 			if(TotalServedMSs-*(CumulativeServingMSs+*(SIndexofUserNumber+i))>=ceil(Ntp*(1-BP))){
// 				printf("%d users will be blocked by Macro cell%d.\n",CumulativeServingMSs[*(SIndexofUserNumber+i)],*(SIndexofUserNumber+i));
// 				for(j=0;j<Ntp;j++)
// 					if(*(temp_dBM+j)==*(SIndexofUserNumber+i)){
// 						*(temp_dBM+j)=-1;
// 						CumulativeServingMSs[*(SIndexofUserNumber+i)]--;
// 						TotalServedMSs--;
// 						*(UnsortedCumulativePower+*(SIndexofUserNumber+i))-=*(P+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(CumulativeDS+*(SIndexofUserNumber+i))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
// 						*(CumulativeBR+*(SIndexofUserNumber+i))-=*(BR+j);
// 						*(UnsortedCumulativeObjective+*(SIndexofUserNumber+i))-=*(SE+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(temp_dBM+j+Ntp)=-1;
// 						BStatus[*(SIndexofUserNumber+i)]=0;
// 					}
// 			}
// 	for(i=0;i<Tbs;i++)
// 		if(*(SIndexofUserNumber+i)>=Tbs-mNbs)
// 			if(TotalServedMSs-*(CumulativeServingMSs+*(SIndexofUserNumber+i))>=ceil(Ntp*(1-BP))){
// 				printf("%d users will be blocked by micro cell%d.\n",CumulativeServingMSs[*(SIndexofUserNumber+i)],*(SIndexofUserNumber+i));
// 				for(j=0;j<Ntp;j++)
// 					if(*(temp_dBM+j)==*(SIndexofUserNumber+i)){
// 						*(temp_dBM+j)=-1;
// 						CumulativeServingMSs[*(SIndexofUserNumber+i)]--;
// 						TotalServedMSs--;
// 						*(UnsortedCumulativePower+*(SIndexofUserNumber+i))-=*(P+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(CumulativeDS+*(SIndexofUserNumber+i))-=*(DS+j*Modulation+*(temp_dBM+j+Ntp));
// 						*(CumulativeBR+*(SIndexofUserNumber+i))-=*(BR+j);
// 						*(UnsortedCumulativeObjective+*(SIndexofUserNumber+i))-=*(SE+(*(SIndexofUserNumber+i)*Ntp*Modulation)+(j*Modulation)+*(temp_dBM+j+Ntp));
// 						*(temp_dBM+j+Ntp)=-1;
// 						BStatus[*(SIndexofUserNumber+i)]=0;
// 					}
// 			}
// 	free(ServingUserNumber);
	for(i=0;i<Tbs;i++)
		printf("%d ", BStatus[i]);
	int* AmountSub=malloc((Tbs-mNbs)*sizeof(int));
	for(i=0;i<Tbs-mNbs;i++){
		if(CumulativeServingMSs[i]==0){
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
			*(UnsortedCumulativePower+i)+=SBP;
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
			*(AmountSub+i)=*(CumulativeDS+i);
			for(j=Tbs-mNbs;j<Tbs;j++)
				if(*(minM+i*Tbs+j)==1)
					*(AmountSub+i)+=*(CumulativeDS+j);
		}
		else{
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
			*(UnsortedCumulativePower+i)+=ABP;
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
			*(AmountSub+i)=*(CumulativeDS+i);
			for(j=Tbs-mNbs;j<Tbs;j++)
				if(*(minM+i*Tbs+j)==1)
					*(AmountSub+i)+=*(CumulativeDS+j);
		}
	}
	for(i=Tbs-mNbs;i<Tbs;i++){
		if(CumulativeServingMSs[i]==0){
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
			*(UnsortedCumulativePower+i)+=mSBP;
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
		}
		else{
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
			*(UnsortedCumulativePower+i)+=mABP;
			printf("BS%d consumes power %lf, status %d, serves %d users\n",i,*(UnsortedCumulativePower+i),BStatus[i],CumulativeServingMSs[i]);
		}
	}
	temp=0;
	for(j=0;j<Tbs;j++)
		if(j<Tbs-mNbs&&(*(UnsortedCumulativePower+j)>(MP+ABP)||(*(AmountSub+j)>DSt))){
			printf("Constraint violation! Program Stops!!\n");
			printf("BS%d Power %lf DS %d\n",j,*(UnsortedCumulativePower+j),*(AmountSub+j)>DSt);
			temp++;
		}
		else if(j>=Tbs-mNbs&&*(UnsortedCumulativePower+j)>(mP+mABP)){ 
			printf("Constraint violation! Program Stops!!\n");
			printf("BS%d Power %lf\n",j,*(UnsortedCumulativePower+j));
			temp++;
		}
	if(temp>0)
		exit(1);
	temp=0;
	for(i=0;i<Tbs;i++)
		temp+=CumulativeServingMSs[i];
	if(TotalServedMSs!=temp){
	printf("6. Inconsistent numbers!%d %d",temp,TotalServedMSs);
	exit(1);
	}
	if(TotalServedMSs!=Ntp-nUser){
		printf("The neccessary number(%d) of total served MSs(%d) isn't achieved!\n",Ntp-nUser,TotalServedMSs);
		exit(1);
	}
	for(i=0;i<Tbs;i++)
		TotalPower+=*(UnsortedCumulativePower+i);
	for(j=0;j<Ntp;j++)
		if(*(temp_dBM+j)>-1)
			NOV+=*(BR+j)/(*(DS+j*Modulation+*(temp_dBM+j+Ntp)));
	double ObjVal=NOV/TotalPower;
	double *BlockedUserSE=malloc((Ntp-TotalServedMSs)*sizeof(double));
	int *BlockedUserIndex=malloc((Ntp-TotalServedMSs)*sizeof(int));
	int *SBlockedUserIndex;
	int *BlockedUsersFirst=malloc(Ntp*sizeof(int));
	l=0;
	for(j=0;j<Ntp;j++)
		if(*(temp_dBM+j)==-1){
			*(BlockedUserSE+l)=0.0;
			*(BlockedUserIndex+l)=-1;
			k=0;
			while((*(AdjacentBSs+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*Ntp+j)==0||BStatus[(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))]==0)&&k<Tbs*Ntp*Modulation-Ntp)
				k+=Ntp;
			if(*(AdjacentBSs+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*Ntp+j)==1&&BStatus[(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))]==1){
				*(BlockedUserSE+l)=*(SE+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*Ntp*Modulation+j*Modulation+(int)(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k)*10-(floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*10)));
				*(BlockedUserIndex+l)=j;
				*(BlockedUsersFirst+j)=(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k));
			}
			l++;
		}
	SBlockedUserIndex=DoubleIndex_q_sort(BlockedUserSE,0,Ntp-TotalServedMSs-1,Ntp-TotalServedMSs,BlockedUserIndex);
	double TObjectiveValue;
	int ReassignmentBS;
	int ReassignmentModulation;
	int CoverMacro;
	for(i=Ntp-TotalServedMSs-1;i>=0;i--){
		if(*(SBlockedUserIndex+i)!=-1){
			TObjectiveValue=ObjVal;
			ReassignmentBS=-1;
			ReassignmentModulation=-1;
			CoverMacro=-1;
			for(j=0;j<Tbs;j++)
				if(*(AdjacentBSs+j*Ntp+*(SBlockedUserIndex+i))==1&&BStatus[j]==1)			
					for(k=0;k<Modulation;k++)
						if((NOV+(*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+k))))/(TotalPower+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k))>=TObjectiveValue){
							if(j<Tbs-mNbs&&*(AmountSub+j)+*(DS+*(SBlockedUserIndex+i)*Modulation+k)<=DSt&&*(UnsortedCumulativePower+j)+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)<=MP+ABP){
								TObjectiveValue=(NOV+(*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+k)))/(TotalPower+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)));
								ReassignmentBS=j;
								ReassignmentModulation=k;
								CoverMacro=-1;
							}
							else if(j>=Tbs-mNbs&&*(UnsortedCumulativePower+j)+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)<=mP+mABP)
								for(l=0;l<Tbs-mNbs;l++)
									if(*(minM+l*Tbs+j)==1&&*(AmountSub+l)+*(DS+*(SBlockedUserIndex+i)+k)<=DSt){
										TObjectiveValue=(NOV+(*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+k)))/(TotalPower+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)));
										ReassignmentBS=j;
										ReassignmentModulation=k;
										CoverMacro=l;
									}
						}
			if(ReassignmentBS!=-1&&ReassignmentModulation!=-1){
				NOV+=*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation));
				TotalPower+=*(P+ReassignmentBS*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
				ObjVal=NOV/TotalPower;
				*(temp_dBM+*(SBlockedUserIndex+i))=ReassignmentBS;
				*(temp_dBM+*(SBlockedUserIndex+i)+Ntp)=ReassignmentModulation;
				CumulativeServingMSs[ReassignmentBS]++;
				TotalServedMSs++;
				*(UnsortedCumulativePower+ReassignmentBS)+=*(P+ReassignmentBS*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
				*(CumulativeDS+ReassignmentBS)+=*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
				*(CumulativeBR+ReassignmentBS)+=*(BR+*(SBlockedUserIndex+i));
				*(UnsortedCumulativeObjective+ReassignmentBS)+=*(SE+ReassignmentBS*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
				if(ReassignmentBS<Tbs-mNbs)
					*(AmountSub+ReassignmentBS)+=*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
				else if(ReassignmentBS>=Tbs-mNbs)
					*(AmountSub+CoverMacro)+=*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
			}
		}
	}
	free(BlockedUserSE);
	free(BlockedUserIndex);
	if(TotalServedMSs<ceil(Ntp*(1-BP))){
		double *BlockedUserSE=malloc((Ntp-TotalServedMSs)*sizeof(double));
		int *BlockedUserIndex=malloc((Ntp-TotalServedMSs)*sizeof(int));
		//int *BlockedUsers=malloc(Ntp*sizeof(int));
		l=0;
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==-1){
				*(BlockedUserSE+l)=0.0;
				*(BlockedUserIndex+l)=-1;
				k=0;
				while((*(AdjacentBSs+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*Ntp+j)==0||BStatus[(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))]==0)&&k<Tbs*Ntp*Modulation-Ntp&&*(BlockedUsersFirst+j)==(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k)))
					k+=Ntp;
				if(*(AdjacentBSs+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*Ntp+j)==1&&BStatus[(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))]==1&&*(BlockedUsersFirst+j)!=(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))){
					*(BlockedUserSE+l)=*(SE+(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*Ntp*Modulation+j*Modulation+(int)(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k)*10-(floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k))*10)));
					*(BlockedUserIndex+l)=j;
					//*(BlockedUsers+j)=(int)floor(*(SIndex+Tbs*Ntp*Modulation-Ntp+j-k));
				}
				l++;
			}
		SBlockedUserIndex=DoubleIndex_q_sort(BlockedUserSE,0,Ntp-TotalServedMSs-1,Ntp-TotalServedMSs,BlockedUserIndex);
		double TObjectiveValue;
		int ReassignmentBS;
		int ReassignmentModulation;
		int CoverMacro;
		for(i=Ntp-TotalServedMSs-1;i>=0;i--){
			if(*(SBlockedUserIndex+i)!=-1){
				TObjectiveValue=ObjVal;
				ReassignmentBS=-1;
				ReassignmentModulation=-1;
				CoverMacro=-1;
				for(j=0;j<Tbs;j++)
					if(*(AdjacentBSs+j*Ntp+*(SBlockedUserIndex+i))==1&&BStatus[j]==1)			
						for(k=0;k<Modulation;k++)
							if((NOV+(*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+k))))/(TotalPower+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k))>=TObjectiveValue){
								if(j<Tbs-mNbs&&*(AmountSub+j)+*(DS+*(SBlockedUserIndex+i)*Modulation+k)<=DSt&&*(UnsortedCumulativePower+j)+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)<=MP+ABP){
									TObjectiveValue=(NOV+(*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+k)))/(TotalPower+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)));
									ReassignmentBS=j;
									ReassignmentModulation=k;
									CoverMacro=-1;
								}
								else if(j>=Tbs-mNbs&&*(UnsortedCumulativePower+j)+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)<=mP+mABP)
									for(l=0;l<Tbs-mNbs;l++)
										if(*(minM+l*Tbs+j)==1&&*(AmountSub+l)+*(DS+*(SBlockedUserIndex+i)+k)<=DSt){
											TObjectiveValue=(NOV+(*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+k)))/(TotalPower+*(P+j*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+k)));
											ReassignmentBS=j;
											ReassignmentModulation=k;
											CoverMacro=l;
										}
							}
				if(ReassignmentBS!=-1&&ReassignmentModulation!=-1){
					NOV+=*(BR+*(SBlockedUserIndex+i))/(*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation));
					TotalPower+=*(P+ReassignmentBS*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
					ObjVal=NOV/TotalPower;
					*(temp_dBM+*(SBlockedUserIndex+i))=ReassignmentBS;
					*(temp_dBM+*(SBlockedUserIndex+i)+Ntp)=ReassignmentModulation;
					CumulativeServingMSs[ReassignmentBS]++;
					TotalServedMSs++;
					*(UnsortedCumulativePower+ReassignmentBS)+=*(P+ReassignmentBS*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
					*(CumulativeDS+ReassignmentBS)+=*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
					*(CumulativeBR+ReassignmentBS)+=*(BR+*(SBlockedUserIndex+i));
					*(UnsortedCumulativeObjective+ReassignmentBS)+=*(SE+ReassignmentBS*Ntp*Modulation+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
					if(ReassignmentBS<Tbs-mNbs)
						*(AmountSub+ReassignmentBS)+=*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
					else if(ReassignmentBS>=Tbs-mNbs)
						*(AmountSub+CoverMacro)+=*(DS+*(SBlockedUserIndex+i)*Modulation+ReassignmentModulation);
				}
			}
		}
		free(BlockedUserSE);
		free(BlockedUserIndex);
		//free(BlockedUsers);
	}
	free(BlockedUsersFirst);
	int *SortedBSforBlockedUser;
	while(TotalServedMSs<ceil(Ntp*(1-BP))){
		int *BSforBlockedUser=malloc(Tbs*sizeof(int));
		for(i=0;i<Tbs;i++){
			*(BSforBlockedUser+i)=0;
			for(j=0;j<nUser;j++)
				if(BStatus[i]==0&&*(AdjacentBSs+i*Ntp+j)==1&&*(temp_dBM+j)==-1)
					*(BSforBlockedUser+i)=*(BSforBlockedUser+i)+1;
		}
		SortedBSforBlockedUser=INT_q_sort(BSforBlockedUser,0,Tbs-1,Tbs);
		for(i=0;i<nUser;i++)
			if(*(AdjacentBSs+*(SortedBSforBlockedUser+Tbs-1)*Ntp+i)==1&&*(temp_dBM+i)==-1){
				k=0;
				double *BlockedSE=malloc(Modulation*sizeof(double));
				int *BlockedSEIndex=malloc(Modulation*sizeof(int));
				for(j=0;j<Modulation;j++){
					*(BlockedSE+j)=*(SE+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+j);
					*(BlockedSEIndex+j)=j;
				}
				DoubleIndex_q_sort(BlockedSE,0,Modulation-1,Modulation,BlockedSEIndex);
				if(*(SortedBSforBlockedUser+Tbs-1)<Tbs-mNbs){
					if(*(AmountSub+*(SortedBSforBlockedUser+Tbs-1))+*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1))<=DSt&&*(UnsortedCumulativePower+*(SortedBSforBlockedUser+Tbs-1))+*(P+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1))<=MP+ABP){
						NOV+=*(BR+i)/(*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1)));
						TotalPower+=*(P+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1));
						if(BStatus[*(SortedBSforBlockedUser+Tbs-1)]==0){
							BStatus[*(SortedBSforBlockedUser+Tbs-1)]=1;
							TotalPower=TotalPower+ABP-SBP;
						}
						ObjVal=NOV/TotalPower;
						*(temp_dBM+i)=*(SortedBSforBlockedUser+Tbs-1);
						*(temp_dBM+i+Ntp)=*(BlockedSEIndex+Modulation-1);
						CumulativeServingMSs[*(SortedBSforBlockedUser+Tbs-1)]++;
						TotalServedMSs++;
						*(UnsortedCumulativePower+*(SortedBSforBlockedUser+Tbs-1))+=*(P+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1));
						*(CumulativeDS+*(SortedBSforBlockedUser+Tbs-1))+=*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1));
						*(CumulativeBR+*(SortedBSforBlockedUser+Tbs-1))+=*(BR+i);
						*(UnsortedCumulativeObjective+*(SortedBSforBlockedUser+Tbs-1))+=*(SE+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1));
						*(AmountSub+*(SortedBSforBlockedUser+Tbs-1))+=*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1));
					}
				}
				else {
					for(j=0;j<Tbs-mNbs;j++)
						if(*(minM+j*Tbs+*(SortedBSforBlockedUser+Tbs-1))==1&&*(AmountSub+j)+*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1))<=DSt&&*(UnsortedCumulativePower)+*(P+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1))<mP+mABP){
							NOV+=*(BR+i)/(*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1)));
							TotalPower+=*(P+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1));
							if(BStatus[*(SortedBSforBlockedUser+Tbs-1)]==0){
								BStatus[*(SortedBSforBlockedUser+Tbs-1)]=1;
								TotalPower=TotalPower-mSBP+mABP;
							}
							ObjVal=NOV/TotalPower;
							*(temp_dBM+i)=*(SortedBSforBlockedUser+Tbs-1);
							*(temp_dBM+i+Ntp)=*(BlockedSEIndex+Modulation-1);
							CumulativeServingMSs[*(SortedBSforBlockedUser+Tbs-1)]++;
							TotalServedMSs++;
							*(UnsortedCumulativePower+*(SortedBSforBlockedUser+Tbs-1))+=*(P+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1));
							*(CumulativeDS+*(SortedBSforBlockedUser+Tbs-1))+=*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1));
							*(CumulativeBR+*(SortedBSforBlockedUser+Tbs-1))+=*(BR+i);
							*(UnsortedCumulativeObjective+*(SortedBSforBlockedUser+Tbs-1))+=*(SE+*(SortedBSforBlockedUser+Tbs-1)*Ntp*Modulation+i*Modulation+*(BlockedSEIndex+Modulation-1));
							*(AmountSub+j)+=*(DS+i*Modulation+*(BlockedSEIndex+Modulation-1));
						}
				}
				free(BlockedSE);
				free(BlockedSEIndex);
			}
		free(BSforBlockedUser);
		SortedBSforBlockedUser=NULL;
	}
	SortedBSforBlockedUser=NULL;
	double *AssignedUserSE=malloc(TotalServedMSs*sizeof(double));
	int *AssignedUserIndex=malloc(TotalServedMSs*sizeof(int));
	int *SortedAssignedUser;
	k=0;
	for(i=0;i<Ntp;i++){
		if(*(temp_dBM+i)>-1){
			*(AssignedUserSE+k)=*(SE+*(temp_dBM+i)*Ntp*Modulation+i*Modulation+*(temp_dBM+i+Ntp));
			*(AssignedUserIndex+k)=i;
			k++;
		}
	}
	SortedAssignedUser=DoubleIndex_q_sort(AssignedUserSE,0,TotalServedMSs-1,TotalServedMSs,AssignedUserIndex);
	for(i=TotalServedMSs-1;i>-1;i--){
		ReassignmentModulation=-1;
		TObjectiveValue=ObjVal;
		CoverMacro=-1;
		for(j=0;j<Modulation;j++){
			if(j!=*(temp_dBM+*(SortedAssignedUser+i)+Ntp)&&(NOV+((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+j)))-((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)))))/(TotalPower+(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+j))-(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))))>TObjectiveValue){
				if(*(temp_dBM+*(SortedAssignedUser+i))<Tbs-mNbs&&*(UnsortedCumulativePower+*(temp_dBM+*(SortedAssignedUser+i)))+*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+j)-*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))<=ABP+MP&&*(AmountSub+*(temp_dBM+*(SortedAssignedUser+i)))+*(DS+*(SortedAssignedUser+i)*Modulation+j)-*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))<=DSt){
					TObjectiveValue=(NOV+((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+j)))-((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)))))/(TotalPower+(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+j))-(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))));
					ReassignmentModulation=j;
					CoverMacro=-1;
				}
				else if(*(temp_dBM+*(SortedAssignedUser+i))>=Tbs-mNbs&&*(UnsortedCumulativePower+*(temp_dBM+*(SortedAssignedUser+i)))+*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+j)-*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))<=mABP+mP){
					for(k=0;k<Tbs-mNbs;k++){
						if(*(minM+k*Tbs+*(temp_dBM+*(SortedAssignedUser+i)))==1&&*(AmountSub+k)+*(DS+*(SortedAssignedUser+i)*Modulation+j)-*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))<=DSt){
							TObjectiveValue=(NOV+((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+j)))-((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)))))/(TotalPower+(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+j))-(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))));
							ReassignmentModulation=j;
							CoverMacro=k;
						}
					}
				}
			}
		}
		if(ReassignmentModulation>-1){
			NOV=NOV+((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+ReassignmentModulation)))-((*(BR+*(SortedAssignedUser+i)))/(*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp))));
			TotalPower=TotalPower+(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+ReassignmentModulation))-(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)));
			ObjVal=NOV/TotalPower;
			*(UnsortedCumulativePower+*(temp_dBM+*(SortedAssignedUser+i)))=*(UnsortedCumulativePower+*(temp_dBM+*(SortedAssignedUser+i)))+(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+ReassignmentModulation))-(*(P+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)));
			*(CumulativeDS+*(temp_dBM+*(SortedAssignedUser+i)))=*(CumulativeDS+*(temp_dBM+*(SortedAssignedUser+i)))+(*(DS+*(SortedAssignedUser+i)*Modulation+ReassignmentModulation))-(*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)));
			*(UnsortedCumulativeObjective+*(temp_dBM+*(SortedAssignedUser+i)))=*(UnsortedCumulativeObjective+*(temp_dBM+*(SortedAssignedUser+i)))+(*(SE+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+ReassignmentModulation))-(*(SE+*(temp_dBM+*(SortedAssignedUser+i))*Ntp*Modulation+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)));
			if(*(temp_dBM+*(SortedAssignedUser+i))<Tbs-mNbs)
				*(AmountSub+*(temp_dBM+*(SortedAssignedUser+i)))=*(AmountSub+*(temp_dBM+*(SortedAssignedUser+i)))+(*(DS+*(SortedAssignedUser+i)*Modulation+ReassignmentModulation))-(*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)));
			else
				*(AmountSub+CoverMacro)=*(AmountSub+CoverMacro)+(*(DS+*(SortedAssignedUser+i)*Modulation+ReassignmentModulation))-(*(DS+*(SortedAssignedUser+i)*Modulation+*(temp_dBM+*(SortedAssignedUser+i)+Ntp)));
			*(temp_dBM+*(SortedAssignedUser+i)+Ntp)=ReassignmentModulation;
		}
	}
	free(AssignedUserSE);
	free(AssignedUserIndex);
	end_tick=clock();
	elapsed=(double)(end_tick-start_tick)/CLOCKS_PER_SEC;	
	int BlockingNumber=0;
	double VerifyingSE=0.0;
	double VerifyingP=0.0;
	temp=0;
	for(i=0;i<Ntp;i++)
		if(*(temp_dBM+i)!=-1&&*(temp_dBM+i+Ntp)!=-1){
			VerifyingSE+=*(BR+i)/(*(DS+i*Modulation+*(temp_dBM+i+Ntp)));
			VerifyingP+=*(P+*(temp_dBM+i)*Ntp*Modulation+i*Modulation+*(temp_dBM+i+Ntp));
			dBM[*(temp_dBM+i)][i][*(temp_dBM+i+Ntp)]=1;
			if(*(AdjacentBSs+*(temp_dBM+i)*Ntp+i)!=1){
				printf("User%d is assigned to infeasible BS%d\n",i,*(temp_dBM+i));
				temp++;
			}
		}
		else if((*(temp_dBM+i)!=-1&&*(temp_dBM+i+Ntp)==-1)||(*(temp_dBM+i)==-1&&*(temp_dBM+i+Ntp)!=-1)){
			printf("Assignment Error!\n");
			exit(1);
		}
		else
			BlockingNumber++;
	
	if(temp>0){
		printf("Program Stops!");
		exit(1);
	}
	temp=0;
	printf("VerifyingP %lf\n",VerifyingP);
	double TEMP=0.0;
	for(i=0;i<Tbs;i++){
		printf("%d ", BStatus[i]);
		TEMP+=*(UnsortedCumulativePower+i);
		temp+=CumulativeServingMSs[i];
		if(BStatus[i]==1&&i<Tbs-mNbs)
			VerifyingP+=ABP;
		else if(BStatus[i]==1&&i>=Tbs-mNbs)
			VerifyingP+=mABP;
		else if(BStatus[i]==0&&i<Tbs-mNbs)
			VerifyingP+=SBP;
		else if(BStatus[i]==0&&i>=Tbs-mNbs)
			VerifyingP+=mSBP;
	}
	printf("VerifyingP %lf\n",VerifyingP);
	if(fabs(VerifyingP-TotalPower)>=0.000001||(VerifyingSE-NOV)>=0.000001||(ObjVal-(VerifyingSE/VerifyingP))>=0.000001){
		printf("\nVerifyingP=%lf TotalPower=%lf VerifyingSE=%lf VerifyingSE/VerifyingP=%lf NOV=%lf ObjVal=%lf TEMP=%lf\n",VerifyingP,TotalPower,VerifyingSE,VerifyingSE/VerifyingP,NOV,ObjVal,TEMP);
		exit(1);
	}
	if(TotalServedMSs!=temp||TotalServedMSs+BlockingNumber!=Ntp){
		printf("6. Inconsistent numbers!%d %d",temp,TotalServedMSs);
		exit(1);
	}
	if(TotalServedMSs<ceil(Ntp*(1-BP))){
		printf("error!!Possible Served MSs:%d. ", TotalServedMSs);
		printf("The required number(%g) of total served MSs(%d) can not be achieved!\n", ceil(Ntp*(1-BP)),Ntp);
		exit(1);
	}
	for(i=0;i<Tbs-mNbs;i++){
		temp=*(CumulativeDS+i);
		for(j=Tbs-mNbs;j<Tbs;j++)
			if(*(minM+i*Tbs+j)==1)
				temp+=*(CumulativeDS+j);
		if(temp!=*(AmountSub+i)){
			printf("Subcarrier numbers don't match! Program Stops!");
			exit(1);
		}
	}
	temp=0;
	for(j=0;j<Tbs;j++)
		if(j<Tbs-mNbs&&(*(UnsortedCumulativePower+j)>(MP+ABP)||(*(AmountSub+j)>DSt))){
			printf("Constraint violation! Program Stops!!\n");
			printf("BS%d Power %lf DS %d\n",j,*(UnsortedCumulativePower+j),*(AmountSub+j)>DSt);
			temp++;
		}
		else if(j>=Tbs-mNbs&&*(UnsortedCumulativePower+j)>(mP+mABP)){ 
			printf("0. Constraint violation! Program Stops!!\n");
			temp++;
			}
	if(temp>0)
		exit(1);
	temp=0;
	for(i=nUser;i<Ntp;i++)
		if(*(temp_dBM+i)==-1){
			printf("User%d must have connection!\n",i);
			temp++;
		}
	if(temp>0){
		printf("Totally %d users have to be connected.\n",temp);
		exit(1);
	}
	free(temp_dBM);
	free(SE);
	free(AmountSub);
	temp=0;
	for(i=0;i<Tbs;i++)
		if(*(CumulativeDS+i)==0&&0-*(CumulativeBR+i)>0.000001){
			printf("BS%d is not supposed to give bandwidth %lf.\n",i,*(CumulativeBR+i));
			temp++;
		}
		else if(*(CumulativeDS+i)!=0&&*(CumulativeBR+i)==0){
			printf("BS%d has no bandwidth to give. DS%d\n",i,*(CumulativeDS+i));
			temp++;
		}
	if(temp>0)
		exit(1);
	int TotalDSs=0;
	double TotalBR=0;
	for(i=0;i<Tbs;i++){
       		printf("\nPower from BS %d: %lf\n",i,*(UnsortedCumulativePower+i));
       		printf("DS from BS %d: %d\n",i,*(CumulativeDS+i));
		printf("BR from BS %d: %lf\n",i,*(CumulativeBR+i));
		printf("Objective from BS %d: %lf\n",i,*(UnsortedCumulativeObjective+i));
       		printf("MSs served by BS %d: %d\n",i,*(CumulativeServingMSs+i));
       		for(j=0;j<Ntp;j++)
       			for(k=0;k<Modulation;k++)
       				if(dBM[i][j][k]==1)
					printf("dBM%d_%d_%d: %d\n",i,j,k,dBM[i][j][k]);
		TotalDSs+=*(CumulativeDS+i);
		TotalBR+=*(CumulativeBR+i);
	}
	double TotalSE=0.0;
	for(i=0;i<Tbs;i++){
		BSSE=0.0;
		BSP=0.0;
       		printf("\nPower from BS %d: %lf\n",i,*(UnsortedCumulativePower+i));
       		for(j=0;j<Ntp;j++)
       			for(k=0;k<Modulation;k++){
       				if(dBM[i][j][k]==1){
       					BSSE+=*(BR+j)/(*(DS+j*Modulation+k));
       					BSP+=*(P+i*Ntp*Modulation+j*Modulation+k);
       				}
       			}
       		printf("Spectrum efficiency from BS %d: %lf\n",i,BSSE);
       		if(*(CumulativeServingMSs+i)>0&&i<Tbs-mNbs)
       			BSP+=ABP;
       		else if(*(CumulativeServingMSs+i)>0&&i>=Tbs-mNbs)
       			BSP+=mABP;
       		else if(*(CumulativeServingMSs+i)==0&&i<Tbs-mNbs)
       			BSP+=SBP;
       		else if(*(CumulativeServingMSs+i)==0&&i>=Tbs-mNbs)
       			BSP+=mSBP;
       		printf("Spectrum-Energy Efficiency of BS %d: %lf\n",i,BSSE/BSP);
       		TotalSE+=BSSE;
	}
  double TranPower = 0.0;
  for(i=0;i<Tbs;i++)
    for(j=0;j<Ntp;j++)
      for(k=0;k<Modulation;k++)
	if(dBM[i][j][k]==1)
	  TranPower += *(P+i*Ntp*Modulation+j*Modulation+k);
  printf("\nTransmit Power: %g", TranPower);
       	printf("\nTotal Power %lf\n",TotalPower);
       	printf("Total spectrum efficiency %lf\n",TotalSE);
       	printf("Total DSs: %d\n",TotalDSs);
	printf("Total BR: %lf\n",TotalBR);
       	printf("Total served MSs: %d\n",TotalServedMSs);
       	printf("New Objective Value: %lf\n",NOV/TotalPower);
       	printf("Running Time: %.9f",elapsed);
	coverage=SCoverageFinder(ptr_dBM,ptrD,Tbs,Ntp,Modulation);
	FILE *HeuristicOut;
	if ((HeuristicOut=fopen("SHeuristicOut.txt", "w")) == NULL)
		printf("\n\nerror!Fail to open file!");
	else
		printf("\n\nOpen SHeuristicOut.txt successfully!\n");
	fprintf(HeuristicOut,"#BS%dMS%dBP%gNU%lf\n",Tbs,Ntp,BP,nUser/Ntp);
	for(i=0;i<Tbs-mNbs;i++)
		fprintf(HeuristicOut," \"BS[%d]\" %d %d %g\n",i,Xbs[i],Ybs[i],*(coverage+i));
	fprintf(HeuristicOut,"\n\n");
	for(i=Tbs-mNbs;i<Tbs;i++)
		fprintf(HeuristicOut," \"BS[%d]\" %d %d %g\n",i,Xbs[i],Ybs[i],*(coverage+i));
	fprintf(HeuristicOut,"\n\n");
        for(i=0;i<Tbs;i++){
		fprintf(HeuristicOut,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
               	for(j=0;j<Ntp;j++)
               		for(k=0;k<Modulation;k++)
				if(dBM[i][j][k]==1){
					fprintf(HeuristicOut,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
					fprintf(HeuristicOut,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW select Mod %d\n",i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],*(DS+j*Modulation+k),*(P+i*Ntp*Modulation+j*Modulation+k),k);
				}
		fprintf(HeuristicOut,"\n\n");
	}
	j=0;
	while(j<Modulation){
		fprintf(HeuristicOut,"\n\n");
		for(i=0;i<Tbs;i++)
			for(k=0;k<Ntp;k++)
				if(dBM[i][k][j]==1)
					fprintf(HeuristicOut,"MS[%d] %d %d DS %d power %g mW select Mod %d\n",k,Xtp[k],Ytp[k],*(DS+k*Modulation+j),*(P+i*Ntp*Modulation+k*Modulation+j),j);
		j++;
	}	
	fprintf(HeuristicOut,"\n\n");
	fprintf(HeuristicOut," plot \"coordinates.txt\" index 0:0 using 2:3 title \"users\" with points, \"SHeuristicOut.txt\" index 0:0 using 2:3:1 notitle with labels, \"SHeuristicOut.txt\" index 0:0 using 2:3:4 title \"Macro Cell\" with circles, \"SHeuristicOut.txt\" index 1:1 using 2:3:1 notitle with labels, \"SHeuristicOut.txt\" index 1:1 using 2:3:4 title \"micro cell\" with circles");
	for(i=0;i<Tbs;i++)
		fprintf(HeuristicOut,", \"SHeuristicOut.txt\" index %d:%d using 5:6 title \"BS[%d]\"",i+2,i+2,i);
	//for(i=0;i<Modulation;i++)
		//fprintf(HeuristicOut,", \"SheuristicOut.txt\" index %d:%d using 2:3 title \"Modulation%d\" with points",Tbs+2+i,Tbs+2+i,i);
	fclose(HeuristicOut);
	//free(coverage);
       	return 0;
}
