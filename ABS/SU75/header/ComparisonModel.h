void cobjective(int Tbs, FILE* PrintToFile,int Ntp, int Modulation, double* BR, int* DS) {
	int count,i,k;
 
//	printf("maximize\n");
	//fprintf(PowerSavingCPLEX,"minimize\n");
	fprintf(PrintToFile,"maximize\n");
// 	for(count=0;count<Ntp;count++)
// 		for(i=0;i<Tbs;i++)
// 			for(k=0;k<Modulation;k++)
// 				if(count==Ntp-1&&i==Tbs-1&&k==Modulation-1)
// 					printf("%lf u%d_%d_%d",BR[count]/DS[count][k],i,count,k);
// 				else 
// 					printf("%lf u%d_%d_%d + ",BR[count]/DS[count][k],i,count,k);
	for(count=0;count<Ntp;count++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(count==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(PrintToFile,"%lf u%d_%d_%d",*(BR+count)/(*(DS+count*Modulation+k)),i,count,k);
				else 
					fprintf(PrintToFile,"%lf u%d_%d_%d + ",*(BR+count)/(*(DS+count*Modulation+k)),i,count,k);
	fprintf(PrintToFile,"\nst\n");
/*	for(count=0; count<Tbs-mNbs; count++){
		for(i=0;i<Ntp;i++)
			for(k=0;k<Modulation;k++)
				if(count==0&&i==0&&k==0)
					fprintf(PowerSavingCPLEX,"%lf dBS + %lf dBM%d_%d_%d + ",SBP*(Tbs-mNbs),P[count][i][k],count,i,k);
				else if(count!=Tbs-mNbs-1&&i==Ntp-1&&k==Modulation-1){
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
					fprintf(PowerSavingCPLEX," + %lf dBA%d + ",ABP-SBP,count);
					//fprintf(PowerSavingCPLEX," + %lf + ",SBP);
				}
				else if(count==Tbs-mNbs-1&&i==Ntp-1&&k==Modulation-1){
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
					fprintf(PowerSavingCPLEX," + %lf dBA%d + ",ABP-SBP,count);						
											
					
					//fprintf(PowerSavingCPLEX," + %lf\n",SBP*Tbs);
				} 
				else
					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[count][i][k],count,i,k);
		
	}*/		
// 	for(count=Tbs-mNbs; count<Tbs; count++)
// 		for(i=0;i<Ntp;i++)
// 			for(k=0;k<Modulation;k++)
// 				if(count==Tbs-mNbs&&i==0&&k==0)
// 					fprintf(PowerSavingCPLEX,"%lf dBS + %lf dBM%d_%d_%d + ",mSBP*mNbs,P[count][i][k],count,i,k);
// 				else if(count!=Tbs-1&&i==Ntp-1&&k==Modulation-1){
// 					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
// 					fprintf(PowerSavingCPLEX," + %lf dBA%d + ",mABP-mSBP,count);
// 					//fprintf(PowerSavingCPLEX," + %lf + ",SBP);
// 				}
// 				else if(count==Tbs-1&&i==Ntp-1&&k==Modulation-1){
// 					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d",P[count][i][k],count,i,k);
// 					fprintf(PowerSavingCPLEX," + %lf dBA%d",mABP-mSBP,count);
// 					//fprintf(PowerSavingCPLEX," + %lf\n",SBP*Tbs);
// 				} 
// 				else
// 					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[count][i][k],count,i,k);
				
	//fprintf(Feasibility, "maximize\n0\nst\n");
return;
}

void ABSconstraint1(int Tbs, int* ptrDS, FILE* PrintToFile, int* minM, int mNbs, int Modulation, int* DS, int DSt, int Ntp) {
	int i,j,k,l,m,n;
	//int o=0;
// 	for(i=0; i<Tbs-mNbs; i++){
// 		//o+=*(Mm+i);
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1){
// 					printf("%d u%d_%d_%d", DS[j][k],i, j,k);
// 					//if(*(Mm+i)>0)
// 					for(l=Tbs-mNbs;l<Tbs;l++)//for(l=Tbs-mNbs+o-*(Mm+i);l<Tbs-mNbs+o;l++)
// 						if(*(minM+i*Tbs+l)==1)
// 							for(m=0;m<Ntp;m++)
// 								for(n=0;n<Modulation;n++)
// 									//if(m==Ntp-1&&n==Modulation-1)
// 										printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);//printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);
// 									//else
// 										//printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);
// 				}		
// 				else 
// 					printf("%d u%d_%d_%d + ",DS[j][k], i, j,k);
// 		}
// 	printf(" - %d t <= 0\n",DSt);	//DSt*200 for 1 second
// 	}
	/*for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(PrintToFile,"%d u%d_%d_%d", *(DS+j*Modulation+k),i, j,k);
					for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									fprintf(PrintToFile," + %d u%d_%d_%d",*(DS+m*Modulation+n),l,m,n);
				}		
				else 
					fprintf(PrintToFile,"%d u%d_%d_%d + ",*(DS+j*Modulation+k), i, j,k);
		}
	fprintf(PrintToFile," - %d t <= 0\n",DSt);	
	}*/
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(PrintToFile,"%d u%d_%d_%d - %d x%d <= 0\n", *(DS+j*Modulation+k),i, j,k, DSt, i);
					for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									if(m==Ntp-1&&n==Modulation-1)
										fprintf(PrintToFile,"%d u%d_%d_%d + %d x%d - %d t <= 0\n",*(DS+m*Modulation+n),l,m,n, DSt, i, DSt);
									else
										fprintf(PrintToFile,"%d u%d_%d_%d + ",*(DS+m*Modulation+n),l,m,n);
				}		
				else 
					fprintf(PrintToFile,"%d u%d_%d_%d + ",*(DS+j*Modulation+k), i, j,k);
		}
		//fprintf(PrintToFile," - %d t <= 0\n",DSt);	
	}
	
	for (i=0; i<Tbs-mNbs; i++){
	  fprintf(PrintToFile, "x%d - t <= 0\n",i);
	  for(l=Tbs-mNbs;l<Tbs;l++)
	    if(*(minM+i*Tbs+l)==1)
	      fprintf(PrintToFile, "x%d + x%d - t = 0\n", i, l);
	}
// 	for(i=0; i<Tbs-mNbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1){
// 					fprintf(PowerSavingCPLEX,"%d dBM%d_%d_%d", DS[j][k],i, j,k);
// 					for(l=Tbs-mNbs;l<Tbs;l++)
// 						if(*(minM+i*Tbs+l)==1)
// 							for(m=0;m<Ntp;m++)
// 								for(n=0;n<Modulation;n++)
// 									fprintf(PowerSavingCPLEX," + %d dBM%d_%d_%d",DS[m][n],l,m,n);
// 				}		
// 				else 
// 					fprintf(PowerSavingCPLEX,"%d dBM%d_%d_%d + ",DS[j][k], i, j,k);
// 		}
// 	fprintf(PowerSavingCPLEX," <= %d\n",DSt);	
// 	}
// 	
// 	for(i=0; i<Tbs-mNbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1){
// 					fprintf(Feasibility,"%d dBM%d_%d_%d", DS[j][k],i, j,k);
// 					for(l=Tbs-mNbs;l<Tbs;l++)
// 						if(*(minM+i*Tbs+l)==1)
// 							for(m=0;m<Ntp;m++)
// 								for(n=0;n<Modulation;n++)
// 									fprintf(Feasibility," + %d dBM%d_%d_%d",DS[m][n],l,m,n);
// 				}		
// 				else 
// 					fprintf(Feasibility,"%d dBM%d_%d_%d + ",DS[j][k], i, j,k);
// 		}
// 	fprintf(Feasibility," <= %d\n",DSt);	
// 	}	
 return;
}//*******************************************************

void cconstraint1(int Tbs, int* ptrDS, FILE* PrintToFile, int* minM, int mNbs, int Modulation, int* DS, int DSt, int Ntp) {
	int i,j,k,l,m,n;
	//int o=0;
// 	for(i=0; i<Tbs-mNbs; i++){
// 		//o+=*(Mm+i);
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1){
// 					printf("%d u%d_%d_%d", DS[j][k],i, j,k);
// 					//if(*(Mm+i)>0)
// 					for(l=Tbs-mNbs;l<Tbs;l++)//for(l=Tbs-mNbs+o-*(Mm+i);l<Tbs-mNbs+o;l++)
// 						if(*(minM+i*Tbs+l)==1)
// 							for(m=0;m<Ntp;m++)
// 								for(n=0;n<Modulation;n++)
// 									//if(m==Ntp-1&&n==Modulation-1)
// 										printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);//printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);
// 									//else
// 										//printf(" + %d u%d_%d_%d",DS[m][n],l,m,n);
// 				}		
// 				else 
// 					printf("%d u%d_%d_%d + ",DS[j][k], i, j,k);
// 		}
// 	printf(" - %d t <= 0\n",DSt);	//DSt*200 for 1 second
// 	}
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(PrintToFile,"%d u%d_%d_%d", *(DS+j*Modulation+k),i, j,k);
					for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									fprintf(PrintToFile," + %d u%d_%d_%d",*(DS+m*Modulation+n),l,m,n);
				}		
				else 
					fprintf(PrintToFile,"%d u%d_%d_%d + ",*(DS+j*Modulation+k), i, j,k);
		}
	fprintf(PrintToFile," - %d t <= 0\n",DSt);	
	}
// 	for(i=0; i<Tbs-mNbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1){
// 					fprintf(PowerSavingCPLEX,"%d dBM%d_%d_%d", DS[j][k],i, j,k);
// 					for(l=Tbs-mNbs;l<Tbs;l++)
// 						if(*(minM+i*Tbs+l)==1)
// 							for(m=0;m<Ntp;m++)
// 								for(n=0;n<Modulation;n++)
// 									fprintf(PowerSavingCPLEX," + %d dBM%d_%d_%d",DS[m][n],l,m,n);
// 				}		
// 				else 
// 					fprintf(PowerSavingCPLEX,"%d dBM%d_%d_%d + ",DS[j][k], i, j,k);
// 		}
// 	fprintf(PowerSavingCPLEX," <= %d\n",DSt);	
// 	}
// 	
// 	for(i=0; i<Tbs-mNbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1){
// 					fprintf(Feasibility,"%d dBM%d_%d_%d", DS[j][k],i, j,k);
// 					for(l=Tbs-mNbs;l<Tbs;l++)
// 						if(*(minM+i*Tbs+l)==1)
// 							for(m=0;m<Ntp;m++)
// 								for(n=0;n<Modulation;n++)
// 									fprintf(Feasibility," + %d dBM%d_%d_%d",DS[m][n],l,m,n);
// 				}		
// 				else 
// 					fprintf(Feasibility,"%d dBM%d_%d_%d + ",DS[j][k], i, j,k);
// 		}
// 	fprintf(Feasibility," <= %d\n",DSt);	
// 	}	
 return;
}//*******************************************************

void cconstraint2(int Tbs,FILE* PrintToFile, int mNbs, int Ntp, int Modulation, double* P, double MP, double mP) {
	int i,j,k;

	/*for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					printf("%lf u%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i, j,k);
				else 
					printf("%lf u%d_%d_%d + ",*(P+i*Ntp*Modulation+j*Modulation+k), i, j,k);
		}
	printf(" - %lf t <= 0\n",MP);
	}
	for(i=Tbs-mNbs; i<Tbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					printf("%lf u%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i, j,k);
				else 
					printf("%lf u%d_%d_%d + ", *(P+i*Ntp*Modulation+j*Modulation+k), i, j,k);
		}
	printf(" - %lf t <= 0\n",mP);
	}*/
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(PrintToFile,"%lf u%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i, j,k);
				else 
					fprintf(PrintToFile,"%lf u%d_%d_%d + ", *(P+i*Ntp*Modulation+j*Modulation+k), i, j,k);
		}
	fprintf(PrintToFile," - %lf t <= 0\n",MP);
	}
	for(i=Tbs-mNbs; i<Tbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(PrintToFile,"%lf u%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i, j,k);
				else 
					fprintf(PrintToFile,"%lf u%d_%d_%d + ",*(P+i*Ntp*Modulation+j*Modulation+k), i, j,k);
		}
	fprintf(PrintToFile," - %lf t <= 0\n",mP);
	}
// 	for(i=0; i<Tbs-mNbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1)
// 					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
// 				else 
// 					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
// 		}
// 	fprintf(PowerSavingCPLEX," <= %lf\n",MP);
// 	}
// 	for(i=Tbs-mNbs; i<Tbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1)
// 					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
// 				else 
// 					fprintf(PowerSavingCPLEX,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
// 		}
// 	fprintf(PowerSavingCPLEX," <= %lf\n",mP);
// 	}
// 	
// 	for(i=0; i<Tbs-mNbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1)
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
// 				else 
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
// 		}
// 		fprintf(Feasibility," <= %lf\n",MP);
// 	}
// 	for(i=Tbs-mNbs; i<Tbs; i++){
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1)
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d", P[i][j][k],i, j,k);
// 				else 
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d + ",P[i][j][k], i, j,k);
// 		}
// 		fprintf(Feasibility," <= %lf\n",mP);
// 	}
 return;
}//*******************************************************

void cconstraint3(int Tbs,FILE* PrintToFile, int nUser, int Modulation) {
 int i,j,k;
// 	 for(i=0; i<nUser; i++)
// 		 for(j=0; j<Tbs; j++) 
// 		 	for(k=0;k<Modulation;k++)
// 				if(j==Tbs-1&&k==Modulation-1)
// 					printf("u%d_%d_%d - t <= 0\n",j,i,k);
// 				else
// 					printf("u%d_%d_%d +",j,i,k);
// 	for(i=0; i<Ntp; i++)
// 		 for(j=0; j<Tbs; j++) 
// 		 	for(k=0;k<Modulation;k++)
// 				if(j==Tbs-1&&k==Modulation-1)
// 					fprintf(NoExistingUser,"u%d_%d_%d - t <= 0\n",j,i,k);
// 				else
// 					fprintf(NoExistingUser,"u%d_%d_%d +",j,i,k);
	for(i=0; i<nUser; i++)
		 for(j=0; j<Tbs; j++) 
		 	for(k=0;k<Modulation;k++)
				if(j==Tbs-1&&k==Modulation-1)
					fprintf(PrintToFile,"u%d_%d_%d - t <= 0\n",j,i,k);
				else
					fprintf(PrintToFile,"u%d_%d_%d +",j,i,k);
	
// 	for(i=0; i<nUser; i++)
// 		 for(j=0; j<Tbs; j++) 
// 		 	for(k=0;k<Modulation;k++)
// 				if(j==Tbs-1&&k==Modulation-1)
// 					fprintf(Feasibility,"dBM%d_%d_%d <= 1\n",j,i,k);
// 				else
// 					fprintf(Feasibility,"dBM%d_%d_%d +",j,i,k);
 return;
}

void cconstraint4(int Tbs, FILE* PrintToFile, int Ntp, int Modulation, double BP) {
 int i,j,k;

//  for(i=0; i<Tbs; i++) 
//   //printf("dbt%d",i);
// 	 for(j=0; j<Ntp; j++)
// 		for(k=0;k<Modulation;k++)
// 			if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
// 				printf("u%d_%d_%d",i,j,k);
// 			else 
// 				printf("u%d_%d_%d +",i,j,k);
// 	 
//   
//  
// printf(" - %.2lf t >= 0\n",Ntp*(1-BP));
	for(i=0; i<Tbs; i++) 
	 	for(j=0; j<Ntp; j++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(PrintToFile,"u%d_%d_%d",i,j,k);
				else 
					fprintf(PrintToFile,"u%d_%d_%d +",i,j,k);
	fprintf(PrintToFile," - %.2lf t >= 0\n",Ntp*(1-BP));
// 	for(i=0; i<Tbs; i++) 
//   //printf("dbt%d",i);
// 	 for(j=0; j<Ntp; j++)
// 		for(k=0;k<Modulation;k++)
// 			if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
// 				fprintf(PowerSavingCPLEX,"dBM%d_%d_%d",i,j,k);
// 			else 
// 				fprintf(PowerSavingCPLEX,"dBM%d_%d_%d +",i,j,k);
// 	 
//   
//  
// fprintf(PowerSavingCPLEX," >= %.2lf\n",Ntp*(1-BP));

// 	for(i=0; i<Tbs; i++) 
// 		for(j=0; j<Ntp; j++)
// 			for(k=0; k<Modulation; k++)
// 				if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
// 					fprintf(Feasibility,"dBM%d_%d_%d",i,j,k);
// 				else 
// 					fprintf(Feasibility,"dBM%d_%d_%d +",i,j,k);
// 	fprintf(Feasibility," >= %.2lf\n",Ntp*(1-BP));
 return;
}

void cconstraint5(int Tbs,FILE* PrintToFile, int Ntp, int Modulation) {
	int i,j,k;

//  	for(i=0;i<Tbs;i++)
// 	 	for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1)
// 					printf("u%d_%d_%d - v%d >= 0\n",i,j,k,i);
// 				else
// 					printf("u%d_%d_%d + ",i,j,k);
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1)
					fprintf(PrintToFile,"u%d_%d_%d - v%d >= 0\n",i,j,k,i);
				else
					fprintf(PrintToFile,"u%d_%d_%d + ",i,j,k);
// 	for(i=0;i<Tbs;i++)
// 	 	for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1)
// 					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d - dBA%d >= 0\n",i,j,k,i);
// 				else
// 					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d + ",i,j,k);
// 				
// 	for(i=0;i<Tbs;i++)
// 	 	for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&k==Modulation-1)
// 					fprintf(Feasibility,"dBM%d_%d_%d - dBA%d >= 0\n",i,j,k,i);
// 				else
// 					fprintf(Feasibility,"dBM%d_%d_%d + ",i,j,k);
return;
}

void cconstraint6(int Tbs,FILE* PrintToFile, int Ntp, int Modulation) {
	int i,j,k;
	int Q=1002;

// 	for(i=0;i<Tbs;i++)
// 	 	for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				//if(j==Ntp-1&&k==Modulation-1)
// 					printf("u%d_%d_%d - v%d <= 0\n",i,j,k,i);
// 				/*else
// 					printf("u%d_%d_%d + ",i,j,k);*/
	for(i=0;i<Tbs;i++)
	 	for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				//if(j==Ntp-1&&k==Modulation-1)
					fprintf(PrintToFile,"u%d_%d_%d - v%d <= 0\n",i,j,k,i);
// 	for(i=0;i<Tbs;i++)
// 	 	for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				fprintf(NoExistingUser,"u%d_%d_%d - v%d <= 0\n",i,j,k,i);
// 			
// 	for(i=0;i<Tbs;i++)
// 	 	for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				fprintf(Feasibility,"dBM%d_%d_%d - dBA%d <= 0\n",i,j,k,i);
 return;
}

void cconstraint7(int Tbs, FILE* PrintToFile, int Ntp, int Modulation) {
	int i,j,k;
	int Q1=1002;

	for(i=0;i<Tbs;i++){
		fprintf(PrintToFile,"v%d - %d b%d <= 0\n",i,Q1,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				fprintf(PrintToFile,"u%d_%d_%d - %d a%d_%d_%d <= 0\n",i,j,k,Q1,i,j,k);
			}
	}
// 	for(i=0;i<Tbs;i++){
// 		printf("v%d - %d b%d <= 0\n",i,Q1,i);
// 		//fprintf(PowerSavingCPLEX,"v%d <= %d t\n",i,Q1);
// 		for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++){
// 				printf("u%d_%d_%d - %d a%d_%d_%d <= 0\n",i,j,k,Q1,i,j,k);
// 				//printf("w%d_%d_%d - %d c%d_%d_%d <= 0\n",i,j,k,Q1,i,j,k);
// 				//fprintf(PowerSavingCPLEX,"u%d_%d <= %d t\n",i,j,Q1);
// 				//fprintf(PowerSavingCPLEX,"w%d_%d <= %d t\n",i,j,Q1);
// 			}
// 	}
return;
}

void cconstraint8(int Tbs, FILE* PrintToFile, int Ntp, int Modulation) {
 int i,j,k;
	int Q2=1002;

// 	for(i=0;i<Tbs;i++){
// 		printf("v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
// 		printf("v%d - t + %d b%d <= %d\n",i,Q2,i,Q2);
// 		//fprintf(PowerSavingCPLEX,"v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
// 		//fprintf(PowerSavingCPLEX,"v%d - t - %d b%d <= %d\n",i,Q2,i,Q2);
// 		for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++){
// 				printf("u%d_%d_%d - t - %d a%d_%d_%d >= -%d\n",i,j,k,Q2,i,j,k,Q2);
// 				printf("u%d_%d_%d - t + %d a%d_%d_%d <= %d\n",i,j,k,Q2,i,j,k,Q2);
// 				//printf("w%d_%d_%d - t - %d c%d_%d_%d >= -%d\n",i,j,k,Q2,i,j,k,Q2);
// 				//printf("w%d_%d_%d - t - %d c%d_%d_%d <= %d\n",i,j,k,Q2,i,j,k,Q2);
// 				//fprintf(PowerSavingCPLEX,"u%d_%d - t - %d a%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
// 				//fprintf(PowerSavingCPLEX,"u%d_%d - t - %d a%d_%d <= %d\n",i,j,Q2,i,j,Q2);
// 				//fprintf(PowerSavingCPLEX,"w%d_%d - t - %d c%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
// 				//fprintf(PowerSavingCPLEX,"w%d_%d - t - %d c%d_%d <= %d\n",i,j,Q2,i,j,Q2);
// 			}
// 	}
	for(i=0;i<Tbs;i++){
		fprintf(PrintToFile,"v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
		fprintf(PrintToFile,"v%d - t + %d b%d <= %d\n",i,Q2,i,Q2);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				fprintf(PrintToFile,"u%d_%d_%d - t - %d a%d_%d_%d >= -%d\n",i,j,k,Q2,i,j,k,Q2);
				fprintf(PrintToFile,"u%d_%d_%d - t + %d a%d_%d_%d <= %d\n",i,j,k,Q2,i,j,k,Q2);
			}
	}
 return;
}

void cconstraint9(int Tbs, FILE* PrintToFile, int mNbs, double ABP, double SBP, int Ntp, int Modulation, double* P, double mABP, double mSBP) {
 int i,j,k;

// 	for(i=0;i<Tbs-mNbs;i++){
// 		printf("%lf v%d + ",ABP-SBP,i);
// 		for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				if(i==Tbs-mNbs-1&&j==Ntp-1&&k==Modulation-1)
// 					printf("%lf u%d_%d_%d",P[i][j][k],i,j,k);
// 				else
// 					printf("%lf u%d_%d_%d + ",P[i][j][k],i,j,k);
// 			
// 	}
// 	for(i=Tbs-mNbs;i<Tbs;i++){
// 		printf(" + %lf v%d",mABP-mSBP,i);
// 		for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				if(i==Tbs-1&&j==Ntp-1&&k==Modulation-1)
// 					printf(" + %lf u%d_%d_%d",P[i][j][k],i,j,k);
// 				else
// 					printf(" + %lf u%d_%d_%d",P[i][j][k],i,j,k);
// 			
// 	}
// 	printf(" + %lf t = 1\n",((Tbs-mNbs)*SBP)+(mNbs*mSBP));
	for(i=0;i<Tbs-mNbs;i++){
		fprintf(PrintToFile,"%lf v%d + ",ABP-SBP,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-mNbs-1&&j==Ntp-1&&k==Modulation-1)
					fprintf(PrintToFile,"%lf u%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i,j,k);
				else
					fprintf(PrintToFile,"%lf u%d_%d_%d + ", *(P+i*Ntp*Modulation+j*Modulation+k),i,j,k);
			
	}
	for(i=Tbs-mNbs;i<Tbs;i++){
		fprintf(PrintToFile," + %lf v%d",mABP-mSBP,i);
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-1&&j==Ntp-1&&k==Modulation-1)
					fprintf(PrintToFile," + %lf u%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i,j,k);
				else
					fprintf(PrintToFile," + %lf u%d_%d_%d",*(P+i*Ntp*Modulation+j*Modulation+k),i,j,k);
			
	}
	fprintf(PrintToFile," + %lf t = 1\n",((Tbs-mNbs)*SBP)+(mNbs*mSBP));

// 	for(j=0;j<Ntp;j++)
// 		for(i=0;i<Tbs;i++)
// 			for(k=0;k<Modulation;k++)
// 				if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d",BR[j]/DS[j][k],i,j,k);
// 				else 
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d + ",BR[j]/DS[j][k],i,j,k);
// 	fprintf(Feasibility," - ");
// 	for(i=0;i<Tbs-mNbs;i++){
// 		fprintf(Feasibility,"%lf dBA%d - ",(ABP-SBP)*limit,i);
// 		for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				if(i==Tbs-mNbs-1&&j==Ntp-1&&k==Modulation-1)
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d",P[i][j][k]*limit,i,j,k);
// 				else
// 					fprintf(Feasibility,"%lf dBM%d_%d_%d - ",P[i][j][k]*limit,i,j,k);
// 	}
// 	for(i=Tbs-mNbs;i<Tbs;i++){
// 		fprintf(Feasibility," - %lf dBA%d",(mABP-mSBP)*limit,i);
// 		for(j=0;j<Ntp;j++)
// 			for(k=0;k<Modulation;k++)
// 				if(i==Tbs-1&&j==Ntp-1&&k==Modulation-1)
// 					fprintf(Feasibility," - %lf dBM%d_%d_%d",P[i][j][k]*limit,i,j,k);
// 				else
// 					fprintf(Feasibility," - %lf dBM%d_%d_%d",P[i][j][k]*limit,i,j,k);
// 			
// 	}
// 	fprintf(Feasibility," - %lf dBS >= 0\n",(((Tbs-mNbs)*SBP)+(mNbs*mSBP))*limit);

 return;
}

void cconstraint10(int Tbs,FILE* PrintToFile, int Ntp, int Modulation, int* AdjacentBSs){
	int i,j,k;

// 	for(i=0;i<Tbs;i++)
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++){
// 				if(k==Modulation-1)
// 					printf("u%d_%d_%d ",i,j,k);
// 				else
// 					printf("u%d_%d_%d + ",i,j,k);
// 			}
// 			printf("- %d t <= 0\n",AdjacentBSs[i][j]);
// 		}
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++){
				if(k==Modulation-1)
					fprintf(PrintToFile,"u%d_%d_%d ",i,j,k);
				else
					fprintf(PrintToFile,"u%d_%d_%d + ",i,j,k);
			}
			fprintf(PrintToFile,"- %d t <= 0\n",*(AdjacentBSs+i*Ntp+j));
		}
// 	for(i=0;i<Tbs;i++)
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++){
// 				if(k==Modulation-1)
// 					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d ",i,j,k);
// 				else
// 					fprintf(PowerSavingCPLEX,"dBM%d_%d_%d + ",i,j,k);
// 			}
// 			fprintf(PowerSavingCPLEX,"<= %d\n",AdjacentBSs[i][j]);
// 		}
		
// 	for(i=0;i<Tbs;i++)
// 		for(j=0;j<Ntp;j++){
// 			for(k=0;k<Modulation;k++){
// 				if(k==Modulation-1)
// 					fprintf(Feasibility,"dBM%d_%d_%d ",i,j,k);
// 				else
// 					fprintf(Feasibility,"dBM%d_%d_%d + ",i,j,k);
// 			}
// 			fprintf(Feasibility,"<= %d\n",AdjacentBSs[i][j]);
// 		}
}

void cconstraint11(int Tbs, FILE* PrintToFile, int nUser, int Ntp, int Modulation){
	int i,j,k;
  if(nUser > Ntp){
    printf("nUser error\n");
    exit(1);
  }
  else if(nUser < Ntp)
	/*for(j=nUser;j<Ntp;j++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(k==Modulation-1&&i==Tbs-1)
					printf("u%d_%d_%d - t = 0\n",i,j,k);
				else
					printf("u%d_%d_%d + ",i,j,k);*/
	for(j=nUser;j<Ntp;j++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(i==Tbs-1&&k==Modulation-1)
					fprintf(PrintToFile,"u%d_%d_%d - t = 0\n",i,j,k);
				else
					fprintf(PrintToFile,"u%d_%d_%d + ",i,j,k);
				
// 	for(j=nUser;j<Ntp;j++)
// 		for(i=0;i<Tbs;i++)
// 			for(k=0;k<Modulation;k++)
// 				if(i==Tbs-1&&k==Modulation-1)
// 					fprintf(Feasibility,"dBM%d_%d_%d = 1\n",i,j,k);
// 				else
// 					fprintf(Feasibility,"dBM%d_%d_%d + ",i,j,k);
}

void ABSbounds(int Tbs, FILE* PrintToFile, int Ntp, int Modulation, int mNbs) {
	int i, j,k;

/*	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				printf("0 <= u%d_%d_%d\n",i,j,k);//dBMi_j
 	for(i=0;i<Tbs;i++) 
		printf("0 <= v%d\n",i);
	//for(i=1;i<=Tbs;i++) 
		printf("0 <= t\n");	*/	
	fprintf(PrintToFile,"bounds\n");
	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				fprintf(PrintToFile,"0 <= u%d_%d_%d\n",i,j,k);
 	for(i=0;i<Tbs;i++) {
		fprintf(PrintToFile,"0 <= v%d\n",i);
		//if (i<Tbs-mNbs)
			fprintf(PrintToFile, "0 <= x%d\n",i);
	}
	fprintf(PrintToFile,"0 <= t\n");

	
// 	fprintf(PowerSavingCPLEX,"dBS=1\n");
// 	
// 	fprintf(Feasibility, "bounds\ndBS=1\n");
	//for(i=1;i<=Tbs;i++) 
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= dBM%d_%d\n",i,j);//dBMi_j
 	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBA%d\n",i);
	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBS%d <= 1\n",i);		
	//for(i=1;i<=Tbs;i++)
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= X%d_%d\n",i,j);
return;
}

void cbounds(int Tbs, FILE* PrintToFile, int Ntp, int Modulation) {
	int i, j,k;

/*	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				printf("0 <= u%d_%d_%d\n",i,j,k);//dBMi_j
 	for(i=0;i<Tbs;i++) 
		printf("0 <= v%d\n",i);
	//for(i=1;i<=Tbs;i++) 
		printf("0 <= t\n");	*/	
	fprintf(PrintToFile,"bounds\n");
	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				fprintf(PrintToFile,"0 <= u%d_%d_%d\n",i,j,k);
 	for(i=0;i<Tbs;i++) 
		fprintf(PrintToFile,"0 <= v%d\n",i);
	fprintf(PrintToFile,"0 <= t\n");

// 	fprintf(PowerSavingCPLEX,"dBS=1\n");
// 	
// 	fprintf(Feasibility, "bounds\ndBS=1\n");
	//for(i=1;i<=Tbs;i++) 
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= dBM%d_%d\n",i,j);//dBMi_j
 	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBA%d\n",i);
	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBS%d <= 1\n",i);		
	//for(i=1;i<=Tbs;i++)
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= X%d_%d\n",i,j);
return;
}

void cspecifyTypes(int Tbs,FILE* PrintToFile, int Ntp, int Modulation) {
	int i,j,k;

 	//printf("binaries\n");
	fprintf(PrintToFile,"binaries\n");
	//fprintf(NoExistingUser,"binaries\n");
/*	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				printf("a%d_%d_%d\n",i,j,k);//dBMi_j
 	for(i=0;i<Tbs;i++) 
		printf("b%d\n",i);*/		
	for(i=0;i<Tbs;i++) 
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				fprintf(PrintToFile,"a%d_%d_%d\n",i,j,k);
 	for(i=0;i<Tbs;i++) 
		fprintf(PrintToFile,"b%d\n",i);
	fprintf(PrintToFile,"end\n");
// 	for(i=0;i<Tbs;i++)
//   		for(j=0;j<Ntp;j++)
//   			for(k=0;k<Modulation;k++)
//    				fprintf(PowerSavingCPLEX,"dBM%d_%d_%d\n",i,j,k);	
// 	for(i=0;i<Tbs;i++)
// 		fprintf(PowerSavingCPLEX,"dBA%d\n",i);
// 	
// 	fprintf(Feasibility, "binaries\n");
// 	for(i=0;i<Tbs;i++)
//   		for(j=0;j<Ntp;j++)
//   			for(k=0;k<Modulation;k++)
//    				fprintf(Feasibility,"dBM%d_%d_%d\n",i,j,k);	
// 	for(i=0;i<Tbs;i++)
// 		fprintf(Feasibility,"dBA%d\n",i);
// 	fprintf(Feasibility, "end\n");
/*	for(i=1;i<=Tbs;i++) 
		fprintf(PowerSavingCPLEX,"dBS%d\n",i);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d\n",i,j);*/

return;
}