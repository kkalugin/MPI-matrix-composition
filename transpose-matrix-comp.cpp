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
	int* ptrA = NULL;
	int* ptrB = NULL;
	int* ptrC = NULL;

    int rank, size;
    int A_size1 = 100, A_size2 = 15, B_size1 = 15, B_size2 = 7;
    int Asize, Bsize, Csize;
    int Bdispls;
    int msgsize;
	int num1, num2;
	int start, end;
	
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

        for(int i = 0; i < B_size2; i++)
            for(int j = 0; j < B_size1; j++)
                B[i * B_size1 + j] = i * j;

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
	
	ptrB = tempB;	
	num1 = (Bdispls + Bsize - 1) / B_size1 - Bdispls  / B_size1;
	start = Bdispls  % B_size1;	
	
	for(int i = 0; i <= num1; i++){
		if(i != num1)
			end = B_size1 - 1;
		else
			end = (Bdispls + Bsize - 1) % B_size1;
		
		ptrC = tempC + Bdispls  / B_size1 + i;
		ptrA = A + start;	
		
		num2 = end - start;		
		for(int j = 0; j < A_size1; j++){						
			for(int k = 0; k <= num2; k++, ptrA++, ptrB++)
				*ptrC +=  *ptrA * *ptrB;			
			ptrB -= num2+1;
			ptrA += A_size2 - num2 - 1;
			ptrC += B_size2;			
		}
		ptrB += num2+1;		
		start = 0;
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
