#define SWAP(x,y) {double t; t = x; x = y; y = t;} 

arrayprint(double *ptr, int BS, int MS){
       int i,j;         
       for(i=0;i<BS;i++){
       //printf("Power Requirement from BS %d:",i);
       for(j=0;j<MS;j++)
       printf("Power Requirement from BS %d to MS %d: %lf \n",i,j,*(ptr+i*(MS)+j));
       //printf("\n");
       //return 0;
       }
}

double * SortedPowerRequirementAssignment(double *PowerRequirment,int Tbs,int Ntp){
       int i,j;
       double *SortedPower = malloc(Tbs*(Ntp) * sizeof(double));
       //printf("%d\n",Tbs*Ntp);
       for(i=0;i<Tbs;i++)
       for(j=0;j<Ntp;j++){
       //*(SortedPower+j*Ntp+i)=*(PowerRequirment+i*Ntp+j);                   
       *(SortedPower+i*Ntp+j)=*(PowerRequirment+(i)*(Ntp)+j);
       //printf("Unsorted Power Requirement from BS %d to MS %d: %lf \n",i+1,j+1,*(SortedPower+i*Ntp+j));
       }
       return SortedPower;                                
}         
                                                  
int * q_sort_2D(double *number, int low, int hight, int MS, int BS){//pass the data array, the starting index and ending index of the data array, and the size of the array
           int *index = malloc(BS*MS * sizeof(int));
           int i,j,k;        
       
//clock_t start_tick, end_tick;           
//double elapsed;
//start_tick = clock();
            
           for(i=0;i<MS;i++)
           for(j=0;j<BS;j++){
           *(index+j*MS+i)=j;
           //printf("index I and J: %d %d ",i,j);
           //printf("content of index MS %d and BS %d %d \n",i,j,*(index+j*MS+i));
           }
           //for(k=0;k<16;k++)
             //printf("%lf ",*(number+k));
           //printf("index size %d\n",BS*(right-left+1));
           for(i=0;i<MS;i++){
           quickSort2D(number,low,hight,MS,i,index,BS);}
           /*for(i=0;i<16;i++){
           printf("%d ",*(index+i));
           printf("%lf ",*(number+i));}*/
           
           
//end_tick = clock();
//elapsed = (double) (end_tick - start_tick) / CLOCKS_PER_SEC; 
//printf("%.100f",elapsed);

           return index;
           }
           
quickSort2D(double *number, int low, int hight, int MS,int MSTerm, int *index, int BS) {  
     //printf("%d %d\n",low,hight);
     if(low < hight) {         
             int i = low;         
             int j = hight + 1;  
             int m;   
             //for(k=0;k<16;k++)
             
             //printf("Start sorting MS %d \n",MSTerm);
             while(1) {       //printf("1. %d %d %d %d\n",i,j,low,hight);                
                      while(i + 1 < BS && *(number+MSTerm+MS*(++i)) < *(number+MSTerm+MS*low)) ;                          
                      while(j -1 > -1 && *(number+MSTerm+MS*(--j)) > *(number+MSTerm+MS*low)) ;              
                      if(i >= j)    {
                           //printf("2. %d %d %d %d\n",low,hight,i,j);             
                           break;            } 
                                   //printf("%lf %lf index %d %d\n",*(number+BS*MS+i),*(number+BS*MS+j),*(index+BS*MS+j),*(index+BS*MS+i));
                                   //printf("SWAP %lf %lf of MS: %d BS: %d and %d\n",*(number+MSTerm+MS*i),*(number+MSTerm+MS*j),MSTerm,*(index+MSTerm+MS*i),*(index+MSTerm+MS*j));
                      SWAP(*(number+MSTerm+MS*i), *(number+MSTerm+MS*j)); 
                      SWAP(*(index+MSTerm+MS*i), *(index+MSTerm+MS*j));
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
             SWAP(*(number+MSTerm+MS*low), *(number+MSTerm+MS*j)); 
             SWAP(*(index+MSTerm+MS*low), *(index+MSTerm+MS*j))
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
             quickSort2D(number, low, j-1, MS, MSTerm, index, BS);  //printf("former ends\n");  printf("later starts\n");        
             quickSort2D(number, j+1, hight, MS, MSTerm, index, BS);/*printf("later ends\n");
             printf("End sorting MS %d \n",MSTerm);       */
             } 
     } 

