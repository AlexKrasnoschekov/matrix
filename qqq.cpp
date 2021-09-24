#include <iostream>
#include "mpi.h"
void Build(int*& lhs, int*& rhs, int n, int m, int k) {
    for (int i = 0; i < n * m; ++i) {
        lhs[i] = i + 1;
    }
    for (int i = 0; i < k * m; ++i) {
        rhs[i] = i - 1;
    }
}
void Print(int*& lhs, int*& rhs, int n, int m, int x) {
	std::cout << "matrix A:" << std::endl;
	for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                        std::cout << lhs[i * m + j] << ' ';
                }
		std::cout << std::endl;
        }
	std::cout << "matrix B:" << std::endl;
	for (int i = 0; i < m; ++i) {
                for (int j = 0; j < x; ++j) {
                        std::cout << rhs[i * x + j] << ' ';
                }
                std::cout << std::endl;
        }
	std::cout << "matrix C:" << std::endl;
}
void GetRhs(int*& matr, int m, int x) {
	int** res = new int*[m];
	for (int i = 0; i < m; ++i) {
		res[i] = new int[x];
	}
	int count = 0;
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < x; ++j) {
			res[i][j] = matr[count];
			count++;
		}
	}
	int** tmp = new int*[x];
        for (int i = 0; i < x; ++i) {
                tmp[i] = new int[m];
        }
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < x; ++j) {
			tmp[j][i] = res[i][j];
		}
	}
	count = 0;
	for (int i = 0; i < x; ++i) {
                for (int j = 0; j < m; ++j) {
                        matr[count] = tmp[i][j];
			count++;
                }
        }
	for (int i = 0; i < m; ++i) {
                delete[] res[i];
        }
	for (int i = 0; i < x; ++i) {
                delete[] tmp[i];
        }
	delete[] tmp;
	delete[] res;
}
int main() {
    MPI_Init(NULL, NULL);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    int n, m, x, *A, *B, *C, *bufA, *bufB, *bufC;
    if (rank == 0) {
        std::cin >> n >> m >> x;
        A = new int[n * m];
        B = new int[x * m];
        Build(A, B, n, m, x);
	Print(A, B, n, m, x);
	GetRhs(B, m, x);
        C = new int[n * x];
        for (int i = 0; i < n * x; ++i) {
            C[i] = 0;
        }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);

    	bufA = new int[n * m / size];
        bufB = new int[x * m / size];
        bufC = new int[n * x / size];
        for (int i = 0; i < n * x / size; ++i) {
                bufC[i] = 0;
        }
    MPI_Scatter(A, n * m / size, MPI_INT, bufA, n * m / size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, x * m / size, MPI_INT, bufB, x * m / size, MPI_INT, 0, MPI_COMM_WORLD);
    for (int q = 0; q < size; ++q) {
        for (int i = 0; i < n / size; ++i) {
            for (int j = 0; j < x / size; ++j) {
                int tmp = 0;
                for (int k = 0; k < m; ++k) {
                    tmp += bufA[i * m + k] * bufB[k + j * m];
            	}
		//std::cout << "rank" << rank << ' ' <<i << ' ' << j << ' ' << tmp << std::endl;
                bufC[(i * x + j + ((size + rank - q) % size) * x / size) % (n * x / size)] = tmp;
	    }
	}
	/*std::cout << rank << ' ' << q << ": ";
	for (int i = 0; i < n * x / size; ++i) {
		std::cout << bufC[i] << ' ';
	}
	std::cout << std::endl;*/
        MPI_Sendrecv_replace(bufB, x * m / size, MPI_INT, (rank + 1) % size, 0, (size + rank - 1) % size, 0, MPI_COMM_WORLD, &status);
    }
    MPI_Gather(bufC, n * x / size, MPI_INT, C, n * x / size, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
	    for (int i = 0; i < n; ++i) {
		    for (int j = 0; j < x; ++j) {
			    std::cout << C[i * x + j] << ' ';
		    }
		    std::cout << std::endl;
	    }
    }
    delete[] bufA;
    delete[] bufB;
    delete[] bufC;
    if (rank == 0) {
        delete[] A;
        delete[] B;
        delete[] C;
    }
    MPI_Finalize();
    return 0;
}
