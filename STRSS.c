#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

//#ifdef DYNAMIC_ALLOC
unsigned char *dna;
//#else
//unsigned char dna[MAX];
//#endif  
unsigned long dna_len;
double total_time;
unsigned long counter;
long int dist1[4][2];
long int dist2[16][2];
long int dist3[64][2];
long int dist4[256][2];
long int dist5[1024][2];
long int dist6[4096][2];
long int repeats1[6];
long int repeats2[6];
long int repeats3[6];
long int repeats4[6];
long int repeats5[6];
long int repeats6[6];
long int repeats[6];
int cut_off; int max_word; int min_word; int p;
unsigned char motif6[7];unsigned char motif5[6];unsigned char motif4[5];
unsigned char motif3[4];unsigned char motif2[3];unsigned char motif1[2];
unsigned char motif[7];
void strs();

unsigned int lookup[20]={0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};
void print_usage() {
    printf("Usage: STRS -l(cut_off)=num -n(minimum word size)=num -m(maximum word size)=num -p(patrial motif)=t/f\n");
}
int main (int argc, char *argv[]) 
{
   // declaring necessary variables to read raw dna file
   FILE *fp1;
   unsigned long i=0; int j=0; 
   double ed = 0.0;
   long seconds, ns; 
   int c;
   opterr=0;
   struct timespec start, finish;
   clock_gettime(CLOCK_REALTIME, &start);
   // opening the data files

   fp1 = fopen(argv[1], "rb");
   if( fp1 == NULL )  {
      perror ("Error opening file");
      return(-1);
   }
   cut_off = 12; min_word =1; max_word =6;
   //Specifying the expected options
   //The three options l, n and m expect numbers as argument
    while ((c = getopt(argc, argv,"l:n:m:p:")) != -1) {
        switch (c) {
             case 'l' : cut_off= atoi(optarg);
                        break;
             case 'n' : min_word = atoi(optarg);
                        break;
             case 'm' : max_word = atoi(optarg); 
                        break;
             case 'p' : p = atoi(optarg); 
                        break;
             default: print_usage(); 
                 exit(EXIT_FAILURE);
        }
    }
    
   printf("CUT_OFF=%d, Min word size=%d, Max word size=%d partial motif =%d \n",cut_off, min_word, max_word,p);
   // checking argument validity

   if (cut_off < 0)
   {
    printf("CUT_OFF value must be positive: Using deafult CUT_OFF value 12\n");
    cut_off = 12;
   }
   else if (min_word < 0)
   {
    printf("Minimum word size must bw within 1 to 6: Using deafult value 1\n");
    min_word = 1;
   }
   else if (max_word < 0)
   {
    printf("Maximum word size must bw within 1 to 6: Using deafult value 6\n");
    max_word = 6;
   }
   else if (min_word > max_word)
   {
    printf("Minimum word size must be less than Maximum word size: Using default value 1 and 6\n");
    min_word = 1;max_word =6;
   }
  if (p != 1 || p !=0)
   printf("Patrial motif can be eiter true(1) or false(0): Using default value =0\n"); 
//#ifdef DYNAMIC_ALLOC
   //dna = (unsigned char *)malloc(3000000000);
//#endif
   //fread(dna, 3000000000 , 1, fp1);
   dna = (unsigned char *)malloc(MAX);
   fread(dna, MAX , 1, fp1);   
   dna_len = strlen(dna); printf("%u",MAX);
   printf("DNA sequene size in bytes:%lu ", dna_len);
   for (i=0;i<6;i++)
   {repeats1[i]=-1;repeats2[i]=-1;repeats3[i]=-1;repeats4[i]=-1;repeats5[i]=-1;repeats6[i]=-1; repeats[i]=-1;}
   
   for (i=0;i<4;i++)
   {dist1[i][0]=-1;dist1[i][1]=-1;}

   for (i=0;i<16;i++)
   {dist2[i][0]=-1;dist2[i][1]=-1;}
   
   for (i=0;i<64;i++)
   {dist3[i][0]=-1;dist3[i][1]=-1;}

   for (i=0;i<256;i++)
   {dist4[i][0]=-1;dist4[i][1]=-1;}

   for (i=0;i<1024;i++)
   {dist5[i][0]=-1;dist5[i][1]=-1;}

   for (i=0;i<4096;i++)
   {dist6[i][0]=-1;dist6[i][1]=-1;}

   clock_gettime(CLOCK_REALTIME, &finish);

   seconds = finish.tv_sec - start.tv_sec;
   ns = finish.tv_nsec - start.tv_nsec;
   total_time = total_time + (double)seconds + (double)ns/(double)100000000;
   printf("\n\nTotal time needed to load sequence and preprocess arguments is : %e \n ", total_time);
 

   clock_gettime(CLOCK_REALTIME, &start);
   
   strs(); // calling short tandem repeat search(strs) 
   clock_gettime(CLOCK_REALTIME, &finish);
   seconds = finish.tv_sec - start.tv_sec;
   ns = finish.tv_nsec - start.tv_nsec;
   total_time = total_time + (double)seconds + (double)ns/(double)100000000;
   printf("\n\nTotal time needed to compute ssr is : %e \n ", total_time);
   
     
   return(0);   
}

unsigned int GetIdxBacktrack_3(unsigned long i)
{
  unsigned long j;
  unsigned long idx = 0;
  idx = 0;
    for(j=i; j<i+3; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
}

unsigned int GetIdxBacktrack_2(unsigned long i)
{
  unsigned long j;
  unsigned long idx = 0;
  for(j=i; j<i+2; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
}

unsigned int GetIdxBacktrack_1(unsigned long i)
{
  unsigned long j;
  unsigned long idx = 0;
  for(j=i; j<i+1; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;    
}


unsigned int GetIdx1(unsigned long i)
{
  unsigned long j;
  static int idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<1; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
  }
   idx = (idx*4+lookup[dna[i+1-1]-'a'])%4;
   return idx;   
}
unsigned int GetIdx2(unsigned long i)
{
  unsigned long j;
  static int idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<2; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
  }
   idx = (idx*4+lookup[dna[i+2-1]-'a'])%16;
   return idx;   
}
unsigned int GetIdx3(unsigned long i)
{
  unsigned long j;
  static int idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<3; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
  }
   idx = (idx*4+lookup[dna[i+3-1]-'a'])%64;
   return idx;   
}
unsigned int GetIdx4(unsigned long i)
{
  unsigned long j;
  static int idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<4; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
  }
   idx = (idx*4+lookup[dna[i+4-1]-'a'])%256;
   return idx;   
}
unsigned int GetIdx5(unsigned long i)
{
  unsigned long j;
  static int idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<5; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
  }
   idx = (idx*4+lookup[dna[i+5-1]-'a'])%1024;
   return idx;   
}
unsigned int GetIdx6(unsigned long i)
{
  unsigned long j;
  static int idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<6; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
  }
   idx = (idx*4+lookup[dna[i+6-1]-'a'])%4096;
   return idx;   
}


void strs()
{
  int j, idx1, idx2, idx3, idx4, idx5, idx6, backtrack=0, forwardtrack=0, d=0, counter=0; int x =0, y=0, end=0;
  long seconds, ns; 
  unsigned long i=0;
  struct timespec start, finish;
  //printf("\nInside strs() : \n dna len = %d\n", dna_len);
  printf("Location:   no. of repeats:    motif:\n" );
  
for(i=0; i<dna_len - 6 + 1; i++) // scanning the sequence one character at a time
{
    // genertaing 4-base codes(index) for k-mers starts at locaion i, for 1<=k<=6.

    idx1 = GetIdx1(i); idx2 = GetIdx2(i); idx3 = GetIdx3(i); 
    idx4 = GetIdx4(i); idx5 = GetIdx5(i); idx6 = GetIdx6(i);
    
    //printf("loc =%d  idx2=%d its prev loc=%d\n",i,idx2, dist2[idx2][0]);
//processing for 6-mer

     if (dist6[idx6][0] == -1) // 6-mer SSR with index = idx6, appearing for first time 
      {
        dist6[idx6][0] = i; // store current location the 6-mer with index = idx6
      } 
     else // not first occurrence, so check for repeats 
      { 
         d = i - dist6[idx6][0];
         if (d < 6) // overlapping occurrence of 6-mer
           {
             dist6[idx6][0] = i; 
             dist6[idx6][1] = -1; // advance the starting point of SSR will deal with atomicity or ; do nothing will report non-atomic SSRs 
           }
         else if (d > 6) // not a consecutive repeating occurrence
           {
              dist6[idx6][0] = i; 
              dist6[idx6][1] = -1; // terminating the current run of 6-mer SSR with index = idx6
           }  
          else // d == 6, consecutive repeating occurrence
           {
               if (dist6[idx6][0] < repeats6[0] + repeats6[5] -6 && idx6 != repeats6[2])// another 6-mer SSR is repeating, or its first time
                  { dist6[idx6][0] = i; dist6[idx6][1] =-1; } // so do not initiate new 6-mer SSR, just advanced its starting
               else 
                  { // initiate new 6-mer SSR 
                    if ( dist6[idx6][1]==-1) // initiation of 6-mer SSR with index = idx6
                    {
                        // check for partial motif 
                        x = repeats6[5] - cut_off; 
                        for (y =0 ; y < x; y++)
                          if ( dna[ repeats6[0] + repeats6[5] - 1 + y] == motif [y]) 
                             repeats6[5] ++;
                        if (repeats6[4] ==1 && repeats6[5] >=cut_off && 6 >= min_word && 6 <= max_word) 
                        // as repeats6 are mutually exclusive, we can print the stored repeats6 
                        {
                           printf("%10ld%10ld%10ld", repeats6[0], repeats6[1], repeats6[3]);
                           printf("%10s\n",motif6);
                           counter++;
                         } 
                         repeats6[0] = dist6[idx6][0];// storing the start location of 6-mer repeats 
                         repeats6[1] = 2; // number of repeats 
                         repeats6[2] = idx6; // storing the 6-mer
                         repeats6[3] = 6; // storing the size of 6-mer, that is 6  
                         repeats6[4] = 1; // setting valid bit
                         repeats6[5] = repeats6[1] * repeats6[3];// computing cut_off 
                         motif6[0] = dna[i]; motif6[1] = dna[i+1];motif6[2] = dna[i+2];motif6[3] = dna[i+3];
                         motif6[4] = dna[i+4]; motif6[5] = dna[i+5]; motif6[6]='\0'; 
                         dist6[idx6][0] = i; 
                         dist6[idx6][1]= 0;
                     }
                    else 
                     {
                        repeats6[1]++; repeats6[5] = repeats6[1] * repeats6[3]; dist6[idx6][0] = i;
                     } 
                  }// end of: initiate and increment 6-mer SSR     
         }// end of: d == 6, consecutive repeating occurrence

       }// end of: not first occurrence, so check for repeats 
   




//processing for 5-mer
 
     if (dist5[idx5][0] == -1) // 5-mer SSR with index = idx5, appearing for first time 
      {
        dist5[idx5][0] = i; // store current location the 5-mer with index = idx5
      } 
     else // not first occurrence, so check for repeats 
      { 
         d = i - dist5[idx5][0];
         if (d < 5) // overlapping occurrence of 5-mer
           {
             dist5[idx5][0] = i; 
             dist5[idx5][1] = -1;;// do nothing 
           }
         else if (d > 5) // not a consecutive repeating occurrence
           {
              dist5[idx5][0] = i; 
              dist5[idx5][1] = -1; // terminating the current run of 5-mer SSR with index = idx5
           }  
          else // d == 5, consecutive repeating occurrence
           {   
               y = i+4; 
               if ((dist5[idx5][0] < repeats5[0] + repeats5[5] -5) && (idx5 != repeats5[2]) || (dist5[idx5][0] >= repeats6[0] && y <= repeats6[0] + repeats6[5] - 1))// another 5-mer SSR is repeating or a 6-mer SSR is completly covering the 5-mer SSR 
                  { dist5[idx5][0] = i; dist5[idx5][1] =-1; } // so do not initiate new 5-mer SSR, just advanced its starting 
               else 
                  { // initiate new 5-mer SSR 
                    if ( dist5[idx5][1]==-1) // initiation of 5-mer SSR with index = idx5
                    {
                        // check for partial motif 
                        end = repeats5[0] + repeats5[5];
                        if (repeats5[4] == 1 && p ==1) 
                        {
                           x = cut_off - repeats5[5]; 
                           for (y = 0 ; y < 5; y++)
                            {
                              //printf("%c %c\n", dna[ end + y ], motif5[y]); 
                              if ( dna[ end + y ] == motif5[y]) 
                                 repeats5[5] ++;
                              else break;  
                            }  
                        } 
                        
                        if (repeats5[4] ==1 && repeats5[5] >=cut_off && 5 >= min_word && 5 <= max_word) 
                         // as repeats5 are mutually exclusive, we can print the stored repeats6 
                        {
                          printf("%10ld%10ld%10ld", repeats5[0], repeats5[1], repeats5[3]);
                          printf("%10s\n",motif5);
                          counter++;
                         } 
                         repeats5[0] = dist5[idx5][0];// storing the start location of 5-mer repeats 
                         repeats5[1] = 2; // number of repeats 
                         repeats5[2] = idx5; // storing the 5-mer
                         repeats5[3] = 5; // storing the size of 5-mer, that is 5  
                         repeats5[4] = 1; // setting valid bit
                         repeats5[5] = repeats5[1] * repeats5[3];// computing cut_off 
                         motif5[0] = dna[i]; motif5[1] = dna[i+1];motif5[2] = dna[i+2];motif5[3] = dna[i+3];
                         motif5[4] = dna[i+4]; motif5[5]='\0'; 
                         dist5[idx5][0] = i; 
                         dist5[idx5][1]= 0;
                     }
                    else 
                     {
                        repeats5[1]++; repeats5[5] = repeats5[1] * repeats5[3]; dist5[idx5][0] = i;
                     } 
                  }// end of: initiate and increment 5-mer SSR     
         }// end of: d == 5, consecutive repeating occurrence

       }// end of: not first occurrence, so check for repeats 





//processing for 4-mer
     //printf("idx4=%d at i=%d\n",idx4,i);
     if (dist4[idx4][0] == -1) // 4-mer SSR with index = idx4, appearing for first time 
      {
        dist4[idx4][0] = i; // store current location the 4-mer with index = idx4
      } 
     else // not first occurrence, so check for repeats 
      { 
         d = i - dist4[idx4][0];
         if (d < 4) // overlapping occurrence of 4-mer
           {
             dist4[idx4][0] = i; 
             dist4[idx4][1] = -1; // terminating the current run of 4-mer SSR with index = idx4;// do nothing 
           }
         else if (d > 4) // not a consecutive repeating occurrence
           {
              dist4[idx4][0] = i; 
              dist4[idx4][1] = -1; // terminating the current run of 4-mer SSR with index = idx4
           }  
          else // d == 4, consecutive repeating occurrence
           {   
               //printf("dist4[idx4][0]=%d and i=%d\n",dist4[idx4][0],i);
               y = i +3;
               if ( ( (dist4[idx4][0] < repeats4[0] + repeats4[5] - 4) && (idx4 != repeats4[2]) ) || (dist4[idx4][0] >= repeats6[0] && y <= repeats6[0] + repeats6[5] - 1) || (dist4[idx4][0] >= repeats5[0] && y <= repeats5[0] + repeats5[5] - 1) )// another 4-mer SSR is repeating or a 6-mer SSR is completly covering the 4-mer SSR or a 5-mer SSR is completly covering the 4-mer SSR 
                  { dist4[idx4][0] = i; dist4[idx4][1] =-1; } // so do not initiate new 4-mer SSR, just advanced its starting 
               else 
                  { // initiate new 4-mer SSR 
                    if ( dist4[idx4][1]==-1) // initiation of 4-mer SSR with index = idx4
                    {
                        if (repeats4[4] ==1 && repeats4[5] >=cut_off && 4 >= min_word && 4 <= max_word) 
                         // as repeats4 are mutually exclusive, we can print the stored repeats4 
                        {
                           printf("%10ld%10ld%10ld", repeats4[0], repeats4[1], repeats4[3]);
                           printf("%10s\n",motif4);
                           counter++;
                         } 
                         repeats4[0] = dist4[idx4][0];// storing the start location of 4-mer repeats 
                         repeats4[1] = 2; // number of repeats 
                         repeats4[2] = idx4; // storing the 4-mer
                         repeats4[3] = 4; // storing the size of 4-mer, that is 4  
                         repeats4[4] = 1; // setting valid bit
                         repeats4[5] = repeats4[1] * repeats4[3];// computing cut_off 
                         motif4[0] = dna[i]; motif4[1] = dna[i+1];motif4[2] = dna[i+2];motif4[3] = dna[i+3];
                         motif4[4] ='\0'; 
                         dist4[idx4][0] = i; 
                         dist4[idx4][1]= 0; //printf("repeats4[0]=%d, i=%d , i+3=%d repeats6[0]=%d repeats6[end]=%d 4-mer motif=%s\n",repeats4[0],i, i+3, repeats6[0], repeats6[0]+ repeats6[5] -1, motif4);
                     }
                    else 
                     {
                        repeats4[1]++; repeats4[5] = repeats4[1] * repeats4[3]; dist4[idx4][0] = i;
                     } 
                  }// end of: initiate and increment 4-mer SSR     
         }// end of: d == 4, consecutive repeating occurrence

       }// end of: not first occurrence, so check for repeats 






//processing for 3-mer
     
     //if ((dna[i] != dna[i+1]) || (dna[i+1] != dna[i+2])) // the 3-mer is non-atomic
     { // process the atomic 3-mer
     if (dist3[idx3][0] == -1) // 3-mer SSR with index = idx3, appearing for first time 
      {
        dist3[idx3][0] = i; // store current location the 3-mer with index = idx3
      } 
     else // not first occurrence, so check for repeats 
      { 
         d = i - dist3[idx3][0];
         if (d < 3) // overlapping occurrence of 3-mer
           {
             dist3[idx3][0] = i; 
             dist3[idx3][1] = -1; // terminating the current run of 3-mer SSR with index = idx3;// do nothing 
           }
         else if (d > 3) // not a consecutive repeating occurrence
           {
              dist3[idx3][0] = i; 
              dist3[idx3][1] = -1; // terminating the current run of 3-mer SSR with index = idx3
           }  
          else // d == 3, consecutive repeating occurrence
           {    
               y = i+2;
               if ( ( (dist3[idx3][0] < repeats3[0] + repeats3[5] - 3) && (idx3 != repeats3[2]) ) || (dist3[idx3][0] >= repeats6[0] && y <= repeats6[0] + repeats6[5] - 1) || (dist3[idx3][0] >= repeats5[0] && y <= repeats5[0] + repeats5[5] - 1) /*|| (dist3[idx3][0] >= repeats4[0] && y <= repeats4[0] + repeats4[5] - 1)*/)// another 3-mer SSR is repeating or a 6-mer SSR is completly covering the 3-mer SSR or a 5-mer SSR is completly covering the 3-mer SSR or a 4-mer SSR is completly covering the 3-mer SSR 
                  { 
                    dist3[idx3][0] = i; dist3[idx3][1] =-1;  // so do not initiate new 3-mer SSR, just advanced its starting 
                  }             
                else 
                  { // initiate new 3-mer SSR 
                    if ( dist3[idx3][1]==-1) // initiation of 3-mer SSR with index = idx3
                    {
                        if (repeats3[4] ==1 && repeats3[5] >=cut_off && 3 >= min_word && 3 <= max_word) 
                         // as repeats3 are mutually exclusive, we can print the stored repeats4 
                        {
                          printf("%10ld%10ld%10ld", repeats3[0], repeats3[1], repeats3[3]);
                          printf("%10s\n",motif3);
                          //printf("%d\n",repeats3[0]);
                          counter++;
                         } 
                         // backtrack the newly founded 3-mer SSR for overlaping 3-mer stared just at previous location 
                         
                         if (dna[i-1] == dna[i-4] && dna[i] == dna[i-3] && dna[i+1] == dna[i-2] ) 
                         {
                             backtrack = GetIdxBacktrack_3(i-4);
                             repeats3[0] = dist3[backtrack][0];// storing the start location of 3-mer repeats 
                             repeats3[1] = 2; // number of repeats 
                             repeats3[2] = backtrack; // storing the 3-mer
                             repeats3[3] = 3; // storing the size of 3-mer, that is 3  
                             repeats3[4] = 1; // setting valid bit
                             repeats3[5] = repeats3[1] * repeats3[3];// computing cut_off 
                             motif3[0] = dna[i-4]; motif3[1] = dna[i-3];motif3[2] = dna[i-2];motif3[3] ='\0'; 
                             dist3[backtrack][0] = i-1; dist3[idx3][0] =i; 
                             dist3[backtrack][1]= 0; dist3[idx3][1] =-1;
                             //printf("backtracking from location %d to %s\n",i,motif3); 
                         }
                         else 
                         {
                             repeats3[0] = dist3[idx3][0];// storing the start location of 3-mer repeats 
                             repeats3[1] = 2; // number of repeats 
                             repeats3[2] = idx3; // storing the 3-mer
                             repeats3[3] = 3; // storing the size of 3-mer, that is 3  
                             repeats3[4] = 1; // setting valid bit
                             repeats3[5] = repeats3[1] * repeats3[3];// computing cut_off 
                             motif3[0] = dna[i]; motif3[1] = dna[i+1];motif3[2] = dna[i+2];motif3[3] ='\0'; 
                             dist3[idx3][0] =i; 
                             dist3[idx3][1] =0;  
                         }
                          
                         //printf("repeats4[0]=%d, i=%d , i+3=%d repeats6[0]=%d repeats6[end]=%d 4-mer motif=%s\n",repeats4[0],i, i+3, repeats6[0], repeats6[0]+ repeats6[5] -1, motif4);
                          
                     }
                    else 
                     {
                        repeats3[1]++; repeats3[5] = repeats3[1] * repeats3[3]; dist3[idx3][0] = i;
                     } 
                  }// end of: initiate and increment 3-mer SSR     
         }// end of: d ==3, consecutive repeating occurrence

       }// end of: not first occurrence, so check for repeats 

} // end: process the atomic 3-mer




//processing for 2-mer
     //printf("idx3=%d at i=%d\n",idx3,i);
     if (dist2[idx2][0] == -1) // 2-mer SSR with index = idx2, appearing for first time 
      {
        dist2[idx2][0] = i; // store current location the 2-mer with index = idx2
      } 
     else // not first occurrence, so check for repeats 
      { 
         d = i - dist2[idx2][0];
         if (d < 2) // overlapping occurrence of 2-mer
           {
             dist2[idx2][0] = i; 
             dist2[idx2][1] = -1; // terminating the current run of 2-mer SSR with index = idx2;// do nothing 
           }
         else if (d > 2) // not a consecutive repeating occurrence
           {
              dist2[idx2][0] = i; 
              dist2[idx2][1] = -1; // terminating the current run of 2-mer SSR with index = idx2
           }  
          else // d == 2, consecutive repeating occurrence
           {   //printf("dist2[idx2][0]=%d and i+1=%d\n",dist2[idx2][0],i+1);
               y = i+1;            
               //printf("dist2[idx2][0]=%d  repeats5[0]%d \n\n",dist2[idx2][0], repeats5[0]);               
               if ( ( (dist2[idx2][0] < repeats2[0] + repeats2[5] - 2) && (idx2 != repeats2[2]) ) || (dist2[idx2][0] >= repeats6[0] && y <= repeats6[0] + repeats6[5] - 1) /*|| ( (dist2[idx2][0] >= repeats5[0]) && (y <= repeats5[0] + repeats5[5] - 1) && (repeats5[5] > cut_off)) || (dist2[idx2][0] >= repeats4[0] && y <= repeats4[0] + repeats4[5] - 1) || (dist2[idx2][0] >= repeats3[0] && y <= repeats3[0] + repeats3[5] - 1)*/)
// another 2-mer SSR is repeating or a 6-mer SSR is completly covering the 2-mer SSR or a 5-mer SSR is completly covering the 2-mer SSR or a 4-mer SSR is completly covering the 2-mer SSR or a 3-mer SSR is completly covering the 2-mer SSR 
                  { dist2[idx2][0] = i; dist2[idx2][1] =-1; } // so do not initiate new 2-mer SSR, just advanced its starting 
               else 
                  { // initiate new 2-mer SSR 
                    if ( dist2[idx2][1]==-1) // initiation of 2-mer SSR with index = idx2
                    {
                       if (repeats2[4] ==1 && repeats2[5] >=cut_off && 2 >= min_word && 2 <= max_word) // as repeats2 are mutually exclusive, we can print the stored repeats4 
                        {
                           printf("%10ld%10ld%10ld", repeats2[0], repeats2[1], repeats2[3]);
                           printf("%10s\n",motif2);
                           counter++;// printf("%d\n",repeats2[0]);
                         } 
                         repeats2[0] = dist2[idx2][0];// storing the start location of 2-mer repeats 
                         repeats2[1] = 2; // number of repeats 
                         repeats2[2] = idx2; // storing the 2-mer
                         repeats2[3] = 2; // storing the size of 2-mer, that is 2  
                         repeats2[4] = 1; // setting valid bit
                         repeats2[5] = repeats2[1] * repeats2[3];// computing cut_off 
                         motif2[0] = dna[i]; motif2[1] = dna[i+1];motif2[2] ='\0'; 
                         dist2[idx2][0] = i; 
                         dist2[idx2][1]= 0; 
                         // backtrack the newly founded 2-mer SSR
                         if (repeats2[4] == 1) backtrack = GetIdxBacktrack_2(repeats2[0] -2);
                         if (backtrack == repeats2[2]) {repeats2[0] = repeats2[0] -2; repeats2[1] ++; repeats2[5] = repeats2[5] + 2; } 
                         // backtrack the newly founded 2-mer SSR for overlaping/cyclic duplicate 2-mer stared just at its previous location 
                         //printf("2-mer ssr found at %d\n",i);
                         if (dna[repeats2[0]-1] == dna[repeats2[0]+1] ) 
                         {
                             backtrack = GetIdxBacktrack_2(repeats2[0]-1);
                             repeats2[0] = dist2[backtrack][0];// storing the start location of 2-mer repeats 
                             repeats2[1] = 2; // number of repeats 
                             repeats2[2] = backtrack; // storing the 2-mer
                             repeats2[3] = 2; // storing the size of 2-mer, that is 2  
                             repeats2[4] = 1; // setting valid bit
                             repeats2[5] = repeats2[1] * repeats2[3];// computing cut_off 
                             motif2[0] = dna[repeats2[0] -1]; motif2[1] = dna[repeats2[0]];motif2[2] ='\0'; 
                             dist2[backtrack][0] = i-1; dist2[idx2][0] =i; 
                             dist2[backtrack][1]= 0; dist2[idx2][1] =-1;
                             //printf("backtracking from location %d to %s at %d\n",i,motif2, repeats2[0]); 
                         }
                         //printf("2-mer SSR starts at %d\n", repeats2[0]); 
                     }// end initiate new 2-mer SSR
                    else 
                     {
                        repeats2[1]++; repeats2[5] = repeats2[1] * repeats2[3]; dist2[idx2][0] = i;
                     } 
                  }// end of: initiate and increment 2-mer SSR     
         }// end of: d ==2, consecutive repeating occurrence

       }// end of: not first occurrence, so check for repeats 





//processing for 1-mer
     //printf("idx1=%d at i=%d\n",idx1,i);
     if (dist1[idx1][0] == -1) // 1-mer SSR with index = idx1, appearing for first time 
      {
        dist1[idx1][0] = i; // store current location the 1-mer with index = idx1
      } 
     else // not first occurrence, so check for repeats 
      { 
         d = i - dist1[idx1][0];
         if (d < 1) // overlapping occurrence of 1-mer
           {
             dist1[idx1][0] = i; 
             dist1[idx1][1] = -1; // terminating the current run of 1-mer SSR with index = idx1;// do nothing 
           }
         else if (d > 1) // not a consecutive repeating occurrence
           {
              dist1[idx1][0] = i; 
              dist1[idx1][1] = -1; // terminating the current run of 1-mer SSR with index = idx1
           }  
          else // d == 1, consecutive repeating occurrence
           {   //printf("dist2[idx2][0]=%d and i+1=%d\n",dist2[idx2][0],i+1);
               y = i;            
               //printf("repeats5[0]=%d  repeats5[0] + repeats5[5] -1 = %d\n",repeats5[0],x);
               
               //if (y<=repeats5[0] + repeats5[5] - 1)  printf("ovelap true i+1=%d x=%d\n\n", i+1,x);
               //else printf("Non-overlap\n"); 
               if ( ( (dist1[idx1][0] < repeats1[0] + repeats1[5] - 1) && (idx1 != repeats1[2]) ) ||  (dist1[idx1][0] >= repeats6[0] && y <= repeats6[0] + repeats6[5] - 1) || (dist1[idx1][0] >= repeats5[0] && y <= repeats5[0] + repeats5[5] - 1) || (dist1[idx1][0] >= repeats4[0] && y <= repeats4[0] + repeats4[5] - 1) || (dist1[idx1][0] >= repeats3[0] && y <= repeats3[0] + repeats3[5] - 1) || (dist1[idx1][0] >= repeats2[0] && y <= repeats2[0] + repeats2[5] - 1))
// another 1-mer SSR is repeating or a 6-mer SSR is completly covering the 1-mer SSR or a 5-mer SSR is completly covering the 1-mer SSR or a 4-mer SSR is completly covering the 1-mer SSR or a 3-mer SSR is completly covering the 1-mer SSR or a 2-mer SSR is completly covering the 1-mer SSR 
                  { dist1[idx1][0] = i; dist1[idx1][1] =-1; } // so do not initiate new 1-mer SSR, just advanced its starting 
               else 
                  { // initiate new 1-mer SSR 
                    if ( dist1[idx1][1]==-1) // initiation of 1-mer SSR with index = idx1
                    {
                       if (repeats1[4] ==1 && repeats1[5] >=cut_off && 1 >= min_word && 1 <= max_word) // as repeats1 are mutually exclusive, we can print the stored repeats4 
                        {
                           printf("%10ld%10ld%10ld", repeats1[0], repeats1[1], repeats1[3]);
                           printf("%10s\n",motif1);
                           counter++;
                         } 
                         
                         repeats1[0] = dist1[idx1][0];// storing the start location of 1-mer repeats 
                         repeats1[1] = 2; // number of repeats 
                         repeats1[2] = idx1; // storing the 1-mer
                         repeats1[3] = 1; // storing the size of 1-mer, that is 1  
                         repeats1[4] = 1; // setting valid bit
                         repeats1[5] = repeats1[1] * repeats1[3];// computing cut_off 
                         motif1[0] = dna[i]; motif1[1] ='\0'; 
                         dist1[idx1][0] = i; 
                         dist1[idx1][1]= 0; 
                         // backtrack the 1-mer SSR: 1 time
                         if (repeats1[4] == 1) backtrack = GetIdxBacktrack_1(repeats1[0] -1);
                         if (backtrack == repeats1[2]) {repeats1[0] = repeats1[0] -1; repeats1[1] ++; repeats1[5] = repeats1[5] + 1; }
                         
                         // backtrack the 1-mer SSR: 2 time
                         if (repeats1[4] == 1) backtrack = GetIdxBacktrack_1(repeats1[0] -1);
                         if (backtrack == repeats1[2]) {repeats1[0] = repeats1[0] -1; repeats1[1] ++; repeats1[5] = repeats1[5] + 1; }

                         // backtrack the 1-mer SSR: 3 time
                         if (repeats1[4] == 1) backtrack = GetIdxBacktrack_1(repeats1[0] -1);
                         if (backtrack == repeats1[2]) {repeats1[0] = repeats1[0] -1; repeats1[1] ++; repeats1[5] = repeats1[5] + 1; }
                         
                         // backtrack the 1-mer SSR: 4 time
                         if (repeats1[4] == 1) backtrack = GetIdxBacktrack_1(repeats1[0] -1);
                         if (backtrack == repeats1[2]) {repeats1[0] = repeats1[0] -1; repeats1[1] ++; repeats1[5] = repeats1[5] + 1; }
//printf("repeats4[0]=%d, i=%d , i+3=%d repeats6[0]=%d repeats6[end]=%d 4-mer motif=%s\n",repeats4[0],i, i+3, repeats6[0], repeats6[0]+ repeats6[5] -1, motif4);
                         //printf("1-mer SSR starts at %d\n", repeats1[0]); 
                     }
                    else 
                     {
                        repeats1[1]++; repeats1[5] = repeats1[1] * repeats1[3]; dist1[idx1][0] = i;
                     } 
                  }// end of: initiate and increment 1-mer SSR     
         }// end of: d ==1, consecutive repeating occurrence

       }// end of: not first occurrence, so check for repeats 

   

}// end scanning 



if (repeats6[4] ==1 && repeats6[5] >=cut_off && 6 >= min_word && 6 <= max_word) // Printing the last stored repeats6 
  {
     printf("%10ld%10ld%10ld", repeats6[0], repeats6[1], repeats6[3]);
     printf("%10s\n",motif6);
     counter++;
   } 

if (repeats5[4] ==1 && repeats5[5] >=cut_off && 5 >= min_word && 5 <= max_word) // Printng the last stored repeats5 
  {
     printf("%10ld%10ld%10ld", repeats5[0], repeats5[1], repeats5[3]);
     printf("%10s\n",motif5);
     counter++;
  } 

if (repeats4[4] ==1 && repeats4[5] >=cut_off && 4 >= min_word && 4 <= max_word) // Printng the last stored repeats4 
  {
     printf("%10ld%10ld%10ld", repeats4[0], repeats4[1], repeats4[3]);
     printf("%10s\n",motif4);
     counter++;
  } 

if (repeats3[4] ==1 && repeats3[5] >=cut_off && 3 >= min_word && 3 <= max_word) // Printng the last stored repeats3 
  {
     printf("%10ld%10ld%10ld", repeats3[0], repeats3[1], repeats3[3]);
     printf("%10s\n",motif3);
     counter++;
  } 


// backtrack the last 2-mer SSR
if (repeats2[4] == 1) {
backtrack = GetIdxBacktrack_2(repeats2[0] -2);
//printf("backtract idx=%d , repeats2[2]=%d\n",backtrack, repeats2[2]);
}
if (backtrack == repeats2[2]) 
{
repeats2[0] = repeats2[0] -2; repeats2[1] ++; 
repeats2[5] = repeats2[5] + 2; //printf("backtracking");
}
if (repeats2[4] ==1 && repeats2[5] >=cut_off && 2 >= min_word && 2 <= max_word) // as repeats2 are mutually exclusive, we can print the stored repeats4 
                        {
                           printf("%10ld%10ld%10ld", repeats2[0], repeats2[1], repeats2[3]);
                           printf("%10s\n",motif2);
                           counter++;
                         } 

if (repeats1[4] ==1 && repeats1[5] >=cut_off && 1 >= min_word && 1 <= max_word) // as repeats1 are mutually exclusive, we can print the stored repeats1 
                        {
                           printf("%10ld%10ld%10ld", repeats1[0], repeats1[1], repeats1[3]);
                           printf("%10s\n",motif1);
                           counter++;
                         } 
printf("Numbers of str are = %d",counter);  
  
}



