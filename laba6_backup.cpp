#include <stdio.h>
#include <mpi.h>



int main(int argc, char *argv[]){
    int* A = NULL;
    int* B = NULL;
    int* C = NULL;
    int* tempC = NULL;
    int* tempB = NULL;
    int* sendcount = NULL;
    int* displs = NULL;

    int rank, size;
    int A_size1 = 100, A_size2 = 15, B_size1 = 15, B_size2 = 7;
    int Asize, Bsize, Csize;
    int Bdispls, Brow, Bcolumn;
    int msgsize;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    msgsize = B_size1 * B_size2 / size ;

    Asize = A_size1 * A_size2;
    Csize = A_size1 * B_size2;

    A = new int [Asize];
    tempC = new int [Csize];
	
    if(rank==0){
        Bsize = B_size1 * B_size2;
        B = new int [Bsize];
        C = new int [Csize];
        sendcount = new int [size];
        displs = new int [size];

        for(int i = 0; i < A_size1; i++)
            for(int j = 0; j < A_size2; j++)
                A[i * A_size2 + j] = i + j;

        for(int i = 0; i < B_size1; i++)
            for(int j = 0; j < B_size2; j++)
                B[i * B_size2 + j] = i * j;

        for(int i = 0; i < Csize; i++){
            tempC[i] = 0;
            C[i] = 0;
        }        

        for(int i = 0; i < size - 1; i++){
            sendcount[i] = msgsize;
            displs[i] = i * msgsize;
        }
        sendcount[size-1] = Bsize - (size - 1) * msgsize;
        displs[size-1] = (size - 1) * msgsize;
	}
	
	
    MPI_Bcast(A, A_size1 * A_size2 , MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(tempC, Csize , MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatter(sendcount, 1, MPI_INT, &Bsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(displs, 1, MPI_INT, &Bdispls, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
    tempB = new int [Bsize];
	
    MPI_Scatterv(B, sendcount, displs, MPI_INT, tempB, Bsize, MPI_INT, 0, MPI_COMM_WORLD);
	
    for(int i = 0; i < Bsize; i++){
        Brow = (Bdispls + i) / B_size2;
        Bcolumn = (Bdispls + i) % B_size2;
        for(int j = 0; j < A_size1; j++)
            tempC[j * B_size2 + Bcolumn] +=  A[j * A_size2 + Brow] * tempB[i];
    }
	
    MPI_Reduce(tempC, C, Csize, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank==0){
        printf("\nMatrix C is:\n");
        for(int i = 0; i < A_size1; i++){
            for(int j = 0; j < B_size2; j++)
                printf("%d\t", C[i * B_size2 + j]);				
            printf("\n");
        }		
        delete []B;
        delete []C;
        delete []sendcount;
        delete []displs;		
    }    
    delete []A;
    delete []tempB;
    delete []tempC;	

    MPI_Finalize();
}
