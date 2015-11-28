#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <mpi.h>
int n, p,z;
int main(int argc, char **argv) {
        int row,col;
        int myrank;

        int *a, *b, *c,*allC, *aRows,*bCols, start, sum, sumdiag;
        int i, j, k;
//my communicators
        MPI_Comm rowcom,colcom;

        n = 4;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD,&p);
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        int myn=n/p;
        row = myrank/p;
        col = myrank-row*p;
        
        a = malloc(myn*sizeof(double));
        b = malloc(myn*sizeof(double));
        c = malloc(myn*sizeof(double));
        aRows = malloc(n*sizeof(double));
        bCols = malloc(n*sizeof(double));

        for(i=0; i<myn; i++) {
                a[i] = 1;
                //printf("a: %i \n",a[i]);
                b[i]=2;
                //printf("b: %i \n",b[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(myrank==0)
                start = MPI_Wtime();

//here I hava all of my parts of my array in the proc


// Creation of my communicators
MPI_Comm_split(MPI_COMM_WORLD, row, myrank,
        &rowcom);
MPI_Comm_split(MPI_COMM_WORLD, col, myrank,
        &colcom);


MPI_Allgather(a, myn, MPI_INT,aRows,n,MPI_INT,rowcom);
MPI_Allgather(b, myn, MPI_INT,bCols,n,MPI_INT,colcom);
MPI_Barrier(MPI_COMM_WORLD);
if(myrank==0){
        for(z=0;z<n;z++)
                printf("arows: %i \n",aRows[z]);
}

        for(j=0; j<2; j++)
                for(i=0; i<2; i++) {

                        sum = 0.;
                        for(k=0; k<n*n; k++)
                                sum += aRows[i*n+k]*bCols[k*n+j];
                        printf("%i",c[i*n+j]);
                        c[i*n+j] = sum;
                }

                MPI_Barrier(MPI_COMM_WORLD);

                if(myrank==0)
                        printf("It took %f seconds to multiply 2 %dx%d matrices.\n",
                                MPI_Wtime()-start, n, n);
                if(myrank==0)
                        allC = malloc(n*n*sizeof(double));
                MPI_Gather(c, row*n, MPI_DOUBLE, allC, row*n, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
                if(myrank==0) {
                        for(i=0, sumdiag=0.; i<n; i++)
                        {
                                sumdiag += allC[i*n+i];
                                printf("vallll :%d \n",allC[i*n+i]);

  }
                        printf("The trace of the resulting matrix is %d\n", sumdiag);
                }

                if(myrank==0)
                        free(allC);
                MPI_Finalize();
                free(a);
                free(b);
                free(c);
        }





