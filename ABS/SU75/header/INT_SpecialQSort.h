#define SWAP(x,y) {double t; t = x; x = y; y = t;} 
#ifndef _2D_QUICKSORT4_H
#define _2D_QUICKSORT4_H


int* Index_q_sort(int number[], int left, int right, int MAX, int *index){//pass the data array, the starting index and ending index of the data array, and the size of the array
           //int *index = malloc((right-left+1) * sizeof(int));
           //int i;
           
//clock_t start_tick, end_tick;           
//double elapsed;
//start_tick = clock();

           //for(i=0;i<MAX;i++)
           //index[i]=i;
           //printf("%d ",index[i]);
           //printf("\n");
           
           Index_quickSort(number,left,right,MAX,index);
           /*for(i=0;i<MAX;i++)
           printf("%d ",number[i]);
           */
           
//end_tick = clock();
//elapsed = (double) (end_tick - start_tick) / CLOCKS_PER_SEC; 
//printf("%.100f",elapsed);

           return index;
           }
           
Index_quickSort(int number[], int left, int right, int MAX,int index[]) {  
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
             Index_quickSort(number, left, j-1,MAX,index);            
             Index_quickSort(number, j+1, right,MAX,index);       
             } 
     } 

#endif
