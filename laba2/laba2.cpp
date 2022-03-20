#include <iostream>
#include <mpi.h>
#include <fstream>
using namespace std;

const int n = 15, m = 10;
double A[n][m], B[m][n], v[n], d[n], C[n][n];

void FillMatrix(double(&AA)[n][m], double(&BB)[m][n]) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			AA[i][j] = rand() % 100 / 2.5;
		}
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			BB[i][j] = rand() % 100 / 2.5;
		}
	}
}
void FillVector(double v1[n]) {
	for (int j = 0; j < n; j++) {
		v1[j] = rand() % 100 / 2.5;
	}
}

void Matrix_Peremnoj() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C[i][j] = 0;
			for (int k = 0; k < m; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}
void Matrix_Peremnoj_na_vector(double A1[n][m], double v1[n]) {
	for (int i = 0; i < n; i++) {
		d[i] = 0;
		for (int j = 0; j < m; j++) {
			d[i] += v1[i] * A1[i][j];
		}
	}
}

void Zapis_v_File() {
	ofstream File1("Matrix_1.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			File1 << A[i][j] << " ";
		}
		File1 << endl;
	}
	File1.close();

	ofstream File2("Matrix_2.txt");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			File2 << B[i][j] << " ";
		}
		File2 << endl;
	}
	File2.close();

	ofstream File4("Vector_1.txt");
	for (int i = 0; i < n; i++)
	{
		File4 << v[i] << endl;
	}
	File4.close();
}

void Zapix_otvetov_v_File() {
	ofstream File3("Matrix_Otvet.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			File3 << C[i][j] << " ";
		}
		File3 << endl;
	}
	File3.close();

	ofstream File5("Vector_Otvet.txt");
	for (int i = 0; i < n; i++)
	{
		File5 << d[i] << endl;
	}
	File5.close();
}

int main()
{
	MPI_Init(NULL, NULL);

	int rank, size;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	FillMatrix(A, B);
	FillVector(v);

	Zapis_v_File();

	Matrix_Peremnoj();

	Matrix_Peremnoj_na_vector(A, v);

	Zapix_otvetov_v_File();





	MPI_Finalize();
}