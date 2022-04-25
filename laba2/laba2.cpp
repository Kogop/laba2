#include <iostream>
#include <mpi.h>
#include <fstream>
using namespace std;

const int root = 0, tag = 0;
const int n = 10, m = 10;
double A[n][m], B[m][n], v[n], d[n], C[n][n], A1[n][m], B1[m][n], v1[n];

void FillMatrix(double(&AA)[n][m], double(&BB)[m][n]) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			AA[i][j] = 1/*rand() % 100 / 2.4*/;
		}
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			BB[i][j] = 1/*rand() % 100 / 2.4*/;
		}
	}
}
void FillVector(double v1[n]) {
	for (int j = 0; j < n; j++) {
		v1[j] = 1/*rand() % 100 / 2.4*/;
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
		for (int j = 0; j < n; j++) {
			cout << CC[i][j] << " ";
		}
		cout << endl;
	}

	//ofstream File5("C:\\Users\\neste\\source\\Repos\\SuperCompModel\\laba2\\laba2\\Vector_Otvet1.txt");
	//for (int i = 0; i < n; i++)
	//{
	//	File5 << d[i] << endl;
	//}
	//File5.close();
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

double k[m];
double l[m];

void vzat_vector_iz_matrix(double(&AA)[n][m], int s1, double(&BB)[m][n], int s2/*, int c*/) {
	// if c == 1 znachit berem stroku if c== 2 znachit berem stolbec
	// s eto nomer stroki ili stolbca kotoriy nuzhno vzat
   // if(c==1){
	for (int j = 0; j < m; j++) {
		k[j] = AA[s1][j];
		//cout << k[j] << endl;
	}
	//  }else if(c==2){

	for (int j = 0; j < m; j++) {
		l[j] = BB[j][s2];
		//cout << l[j] << endl;
	}
	//}
}

double peremnoj_vector_na_vector(double(&kk)[m], double(&ll)[m]) {
	double h = 0;
	for (int i = 0; i < m; i++) {
		h += kk[i] * ll[i];
	}
	return h;
}

double Temp[] = { 0,0,0 };

int main() {
	MPI_Init(NULL, NULL);

	int rank, size, limit, end, end_1_otprav = 0, end_1_priem = 0;
	end = 0;

	MPI_Status status;
	//MPI_Request request;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int tag1 = 1, tag2 = 2, tag3 = 0;
	int isDone = 0;

	if (n < size) {
		limit = n + 1;  // if kolichestvo processov bolwe razmera matrix
	}
	else { limit = size; }  // if kolichestvo processov menwe ili ravno razmeru matrix

	if (rank == 0) {
		FillMatrix(A, B);
		//FillVector(v);
		Zapis_v_File();
		read_Vector();
		read_Matrix();

		for (int end = 0; end < n * n; end += end_1_otprav)  // 4et  poka viglyadit nepravilno
		{
			end_1_otprav = 0;
			cout << "some part of all ranks" << " started with " << end << endl;
			for (int i = 0; i < limit - 1; i++)
			{
				for (int j = 0; j < limit - 1; j++)
				{
					int ii = end + i;
					int jj = end + j;
					vzat_vector_iz_matrix(A1, ii, B1, jj/*,1*/);
					//tag1 = i;
					MPI_Send(k, m, MPI_DOUBLE, j + 1, tag1, MPI_COMM_WORLD);
					// vzat_vector_iz_matrix(B1,j/*,2*/);
					 //tag2 = j;
					MPI_Send(l, m, MPI_DOUBLE, j + 1, tag2, MPI_COMM_WORLD);
					cout << "Otpr 1 rannk = " << rank << " i = " << i << " j = " << j << endl;
					fflush(stdout);
					end_1_otprav++;
				}
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				// tag3 = j;
				MPI_Recv(&(Temp[0]), 3, MPI_DOUBLE, MPI_ANY_SOURCE, tag3, MPI_COMM_WORLD, &status);
				int h = Temp[1] - 1;
				int g = Temp[2];
				C[h][g] = Temp[0];
				cout << "Priem 2 rank = " << j + 1 << " C = " << C[h][g] << endl;
			}
		}
		
		//MPI_Barrier(MPI_COMM_WORLD);
		Zapix_otvetov_v_File(C/*,d*/);

	}
	else if (rank < limit) {
		
		for (int end = 0; end < n*n; end += end_1_priem)  // 4et  poka viglyadit nepravilno
		{
			end_1_priem = 0;
			for (int j = 0; j < limit - 1; j++) {
				//tag3 = 0;
				MPI_Recv(&(k[0]), m, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);
				MPI_Recv(&(l[0]), m, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);
				cout << "Priem 1 rank = " << rank << " j = " << j << endl;
				//Matrix_Peremnoj(A1, B1);
			   // C[i][j] = peremnoj_vector_na_vector(k,l);
				Temp[0] = peremnoj_vector_na_vector(k, l);
				Temp[1] = rank;
				Temp[2] = j;
				MPI_Send(&Temp, 3, MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD);
				end_1_priem++;
			}
		}
		
		cout << rank << " all counting is done. Writing the answers" << endl;
	}
	//   MPI_Barrier(MPI_COMM_WORLD);
	 // if(rank == 1){
	  // MPI_Send(C, n*n, MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD);
	  // }
   // uznavat' v kakom range zapisalsya element i iz nego uze poluchat vse

	MPI_Finalize();
	return 1;
}
