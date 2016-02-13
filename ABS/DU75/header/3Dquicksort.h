#define SWAP(x,y) {double t; t = x; x = y; y = t;} 
#ifndef _3D_QUICKSORT_H
#define _3D_QUICKSORT_H

//#include "2Dquicksort4.h"
//quickSort2D(double *number, int low, int hight, int MS, int MSTerm, int *index, int BS)
double * q_sort_3D(double *Power, int Tbs, int Ntp, int Modulation){
	int i,j,k,m,n;
	double *BSMod=malloc(Ntp*(Tbs*Modulation)*sizeof(double));
	double *IndexBSMod=malloc(Ntp*(Tbs*Modulation)*sizeof(double));
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++){
				*(IndexBSMod+j+i*(Ntp*Modulation)+k*Ntp)=i+(k*0.1);
				//printf("%lf\n",*(IndexBSMod+j*(Tbs*Modulation)+i*Modulation+k));
				*(BSMod+i*(Ntp*Modulation)+k*Ntp+j)=*(Power+i*Ntp*Modulation+j*Modulation+k);
				//printf("%lf\n",*(BSMod+j*(Tbs*Modulation)+i*Modulation+k));
				//printf("%d MS%d and BS %d use modulation %d transmitting power %lf. \n",i*Ntp*Modulation+j*Modulation+k,j,floor(*(IndexBSMod+j*(Tbs*Modulation)+i*Modulation+k)),*(IndexBSMod+j*(Tbs*Modulation)+i*Modulation+k)*10-floor(*(IndexBSMod+j*(Tbs*Modulation)+i*Modulation+k))*10,*(BSMod+j*(Tbs*Modulation)+i*Modulation+k));
			}
	/*for(n=0;n<Ntp;n++){
		for(m=0;m<Tbs*Modulation;m++){
			//printf("%d MS%d and BS%d use modulation %d transmitting power %lf. \n",n*Tbs*Modulation+m,n,(int)floor(*(IndexBSMod+n*(Tbs*Modulation)+m)),(int)(*(IndexBSMod+n*(Tbs*Modulation)+m)*10-floor(*(IndexBSMod+n*(Tbs*Modulation)+m))*10),*(BSMod+n*(Tbs*Modulation)+m));
			printf("%d MS%d %lf transmitting power %lf. \n",n*Tbs*Modulation+m,n,(*(IndexBSMod+m*Ntp+n)),*(BSMod+m*Ntp+n));
		}
		printf("\n");
	}*/
	//for(m=0;m<Tbs*Modulation;m++)
	//	printf("%d MS%d %lf transmitting power %lf. \n",m,0,(*(IndexBSMod+m*Ntp+0)),*(BSMod+m*Ntp+0));
	printf("\n\n");
	for(i=0;i<Ntp;i++)
		DoubleQuickSort2D(BSMod,0,Tbs*Modulation-1,Ntp,i,IndexBSMod,Tbs*Modulation);
	//for(m=0;m<Tbs*Modulation;m++)
	//	printf("%d MS%d %lf transmitting power %lf. \n",m,0,(*(IndexBSMod+m*Ntp+0)),*(BSMod+m*Ntp+0));
	/*for(n=0;n<Ntp;n++){
		for(m=0;m<Tbs*Modulation;m++){
			printf("%d MS%d and BS%d use modulation %d transmitting power %lf. \n",n*Tbs*Modulation+m,n,(int)floor(*(IndexBSMod+m*Ntp+n)),(int)(*(IndexBSMod+m*Ntp+n)*10-floor(*(IndexBSMod+m*Ntp+n))*10),*(BSMod+m*Ntp+n));
			//printf("%d MS%d %lf transmitting power %lf. \n",n*Tbs*Modulation+m,n,(*(IndexBSMod+n*(Tbs*Modulation)+m)),*(BSMod+n*(Tbs*Modulation)+m));
		}
		printf("\n");
	}*/
return IndexBSMod;
}
DoubleQuickSort2D(double *PArray, int lowest, int hightest, int SS,int SSTerm, double *PAIndex, int Length) {  
     //printf("%d %d\n",low,hight);
     if(lowest < hightest) {         
             int i = lowest;         
             int j = hightest + 1;  
             int m;   
             //for(k=0;k<16;k++)
             
             //printf("Start sorting MS %d \n",MSTerm);
             while(1) {       //printf("1. %d %d %d %d\n",i,j,low,hight);                
                      while(i + 1 < Length && *(PArray+SSTerm+SS*(++i)) < *(PArray+SSTerm+SS*lowest)) ;                          
                      while(j -1 > -1 && *(PArray+SSTerm+SS*(--j)) > *(PArray+SSTerm+SS*lowest)) ;              
                      if(i >= j)    {
                           //printf("2. %d %d %d %d\n",low,hight,i,j);             
                           break;            } 
                                   //printf("%lf %lf index %d %d\n",*(number+BS*MS+i),*(number+BS*MS+j),*(index+BS*MS+j),*(index+BS*MS+i));
                                   //printf("SWAP %lf %lf of MS: %d BS: %lf and %lf\n",*(PArray+SSTerm+SS*i),*(PArray+SSTerm+SS*j),SSTerm,*(PAIndex+SSTerm+SS*i),*(PAIndex+SSTerm+SS*j));
                      SWAP(*(PArray+SSTerm+SS*i), *(PArray+SSTerm+SS*j)); 
                      SWAP(*(PAIndex+SSTerm+SS*i), *(PAIndex+SSTerm+SS*j));
                                   //printf("After SWAP %lf %lf of MS: %d BS: %d and %d\n",*(number+MSTerm+MS*i),*(number+MSTerm+MS*j),MSTerm,*(index+MSTerm+MS*i),*(index+MSTerm+MS*j));
                      //index[i]=j;
                      //index[j]=i;
                                   /*printf("inner new array:");
                                   for(m=0;m<16;m++)
                                   printf("%lf ",*(number+m));        
                                   printf("\n");
                                   printf("inner new index array:");
                                   for(m=0;m<16;m++)
                                   printf("%d ",*(index+m));        
                                   printf("\n");*/
                                   
                      }         
             SWAP(*(PArray+SSTerm+SS*lowest), *(PArray+SSTerm+SS*j)); 
             SWAP(*(PAIndex+SSTerm+SS*lowest), *(PAIndex+SSTerm+SS*j))
             //index[left]=j;
             //index[j]=left;
                             /*printf("After SWAP %lf %lf index %d %d\n",*(number+MSTerm+MS*low),*(number+MSTerm+MS*j),*(index+MSTerm+MS*j),*(index+MSTerm+MS*low));
                             printf("outer new array:");
                             for(m=0;m<16;m++)
                             printf("%lf ",*(number+m));        
                             printf("\n");
                             printf("outer new index array:");
                             for(m=0;m<16;m++)
                             printf("%d ",*(index+m));        
                             printf("\n");
                             printf("former starts\n");
                             printf("3. %d %d %d %d\n",low,hight,i,j );*/
             DoubleQuickSort2D(PArray, lowest, j-1, SS, SSTerm, PAIndex, Length);  //printf("former ends\n");  printf("later starts\n");        
             DoubleQuickSort2D(PArray, j+1, hightest, SS, SSTerm, PAIndex, Length);/*printf("later ends\n");
             printf("End sorting MS %d \n",MSTerm);       */
             } 
     } 
QuickSort(double *number, int left, int right, int MAX,double *index) {  
     if(left < right) {         
             int i = left;         
             int j = right + 1;  
             int m;       
             while(1) {                       
                      while(i + 1 < MAX && number[++i] < number[left]) ;                          
                      while(j -1 > -1 && number[--j] > number[left]) ;              
                      if(i >= j)                 
                           break;             
                      SWAP(number[i], number[j]); 
                      SWAP(index[i],index[j]);
                      /*printf("SWAP %d %d index %d %d\n",number[i],number[j],j,i);
                      index[i]=j;
                      index[j]=i;
                      printf("inner new array:");
                      for(m=0;m<MAX;m++)
                      printf("%d ",number[m]);        
                      printf("\n");
                      printf("inner new index array:");
                      for(m=0;m<MAX;m++)
                      printf("%d ",index[m]);        
                      printf("\n");*/
                      }         
             SWAP(number[left], number[j]); 
             SWAP(index[left],index[j])
             //index[left]=j;
             //index[j]=left;
             /*printf("SWAP %d %d index %d %d\n",number[left],number[j],index[j],index[left]);
              printf("outer new array:");
                       for(m=0;m<MAX;m++)
                       printf("%d ",number[m]);        
                       printf("\n");
                     printf("outer new index array:");
                      for(m=0;m<MAX;m++)
                      printf("%d ",index[m]);        
                      printf("\n");*/
             QuickSort(number, left, j-1,MAX,index);            
             QuickSort(number, j+1, right,MAX,index);       
             } 
     } 
#endif
