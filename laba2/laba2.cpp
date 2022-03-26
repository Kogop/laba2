﻿#include <iostream>
#include <mpi.h>
#include <fstream>
using namespace std;

const int root = 0, tag = 0;
const int n = 15, m = 10;
double A[n][m], B[m][n], v[n], d[n], C[n][n], A1[n][m], B1[m][n], v1[n];

void FillMatrix(double(&AA)[n][m], double(&BB)[m][n]) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			AA[i][j] = rand() % 100 / 2.4;
		}
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			BB[i][j] = rand() % 100 / 2.4;
		}
	}
}
void FillVector(double v1[n]) {
	for (int j = 0; j < n; j++) {
		v1[j] = rand() % 100 /2.4;
	}
}

void Matrix_Peremnoj(double(&AA)[n][m], double(&BB)[m][n]) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C[i][j] = 0;
			for (int k = 0; k < m; k++) {
				C[i][j] += AA[i][k] * BB[k][j];
			}
		}
	}
}
void Matrix_Peremnoj_na_vector(double(&AA)[n][m], double(&vv)[n]) {
	for (int i = 0; i < n; i++) {
		d[i] = 0;
		for (int j = 0; j < m; j++) {
			d[i] += vv[i] * AA[i][j];
		}
	}
}

void Zapis_v_File() {
	ofstream File1("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Matrix_1.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			File1 << A[i][j] << " ";
		}
		File1 << endl;
	}
	File1.close();

	ofstream File2("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Matrix_2.txt");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			File2 << B[i][j] << " ";
		}
		File2 << endl;
	}
	File2.close();

	ofstream File4("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Vector_1.txt");
	for (int i = 0; i < n; i++)
	{
		File4 << v[i] << endl;
	}
	File4.close();
}

void Zapix_otvetov_v_File(double(&CC)[n][n]) {
	ofstream File3("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Matrix_Otvet1.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			File3 << CC[i][j] << " ";
		}
		File3 << "\n";
	}
	File3.close();
	cout << " c zapisannaya --" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			cout << CC[i][j] << " ";
		}
		cout << endl;
	}

	ofstream File5("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Vector_Otvet1.txt");
	for (int i = 0; i < n; i++)
	{
		File5 << d[i] << endl;
	}
	File5.close();
}

void read_Vector() {
	ifstream File5("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Vector_1.txt");
	for (int i = 0; i < n; i++) {
		File5 >> v1[i];
		//  cout << DD[i]<<endl;
	}
	File5.close();
}

void read_Matrix() {
	ifstream File1("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Matrix_1.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			File1 >> A1[i][j];
			//cout << AA[i][j] << " ";
		}
		// cout << endl;
	}
	File1.close();

	ifstream File2("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Matrix_2.txt");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			File2 >> B1[i][j];
			//cout << AA[i][j] << " ";
		}
		// cout << endl;
	}
	File2.close();
}

int main()
{
	MPI_Init(NULL, NULL);

	int rank, size;
	MPI_Status status;
	//MPI_Request request;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int tag1 = 1, tag2 = 2, tag3 = 3;
	int isDone = 0;

	if (rank == 0) {
		FillMatrix(A, B);
		FillVector(v);
		Zapis_v_File();
		read_Vector();
		read_Matrix();

		MPI_Send(A1, n*m, MPI_DOUBLE, 1, tag1, MPI_COMM_WORLD);

		MPI_Send(B1, m*n, MPI_DOUBLE, 1, tag1, MPI_COMM_WORLD);

		MPI_Send(v1, 1, MPI_DOUBLE, 2, tag2, MPI_COMM_WORLD);
		MPI_Send(A1, 1, MPI_DOUBLE, 2, tag2, MPI_COMM_WORLD);

		//isDone++;
		cout << rank << "th rank is done it's work flawlessly" << endl;
	}
	else if (rank == 1) {
		MPI_Recv(&(A1[0][0]), n*m, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);

		MPI_Recv(&(B1[0][0]), m*n, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);

		Matrix_Peremnoj(A1, B1);
		MPI_Send(C, n*n, MPI_DOUBLE, 3, tag3, MPI_COMM_WORLD);
		//isDone++;
		cout << rank << "st rank is done it's work flawlessly" << endl;
	}
	else if (rank == 2)
	{
		MPI_Recv(&(v1[0]), n, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);
		MPI_Recv(&(A1[0][0]), n*m, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);

		Matrix_Peremnoj_na_vector(A1, v1);

		MPI_Send(d, n, MPI_DOUBLE, 3, tag3, MPI_COMM_WORLD);
		//isDone++;
		cout << rank << "nd rank is done it's work flawlessly" << endl;
	}
	else if (rank == 3)
	{
		MPI_Recv(&(C[0][0]), n*n, MPI_DOUBLE, 1, tag3, MPI_COMM_WORLD, &status);
		MPI_Recv(&(d[0]), n, MPI_DOUBLE, 2, tag3, MPI_COMM_WORLD, &status);
		//cout << " c poluchennaya --" << endl;
		//for (int i = 0; i < n; i++)
		//{
		//	for (int j = 0; j < m; j++) {
		//		cout << C[i][j] << " ";
		//	}
		//	cout << endl;
		//}

		Zapix_otvetov_v_File(C);
		//isDone++;
		cout << rank << "rd rank is done it's work flawlessly" << endl;
	}
	else if(rank!=0)
	{
		for (int rank = 0; rank < 3; rank++)
		{

		}
	}
	//cout << isDone << endl;

	MPI_Finalize();
}