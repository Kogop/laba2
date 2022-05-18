#include <iostream>
#include <mpi.h>
#include <fstream>
//#include <math.h>
#include <cmath>
#include <stdlib.h>
using namespace std;

const int root = 0, tag = 0;
const int n = 100, m = 150;
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
/*
void FillVector(double v1[n]) {
	for (int j = 0; j < n; j++) {
		v1[j] = 1/*rand() % 100 / 2.4;
	}
}
*/

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
			File3 << CC[i][j] << "\t";
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
int Mesto[] = { 0,0 };
int N = 0, M = 0, Nn = 0, Mm = 0, limit_1 = 0;
int kol_strok_v_posled_str_bloke = 0;
int kol_stolb_v_posled_stolb_bloke = 0;
bool is_n_menwe_size;

int main() {
	MPI_Init(NULL, NULL);
	double starttime, endtime;
	starttime = MPI_Wtime();

	int rank, size, limit, end, end_1_otprav = 0, end_1_priem = 0, h = 0, g = 0;
	end = 0;

	MPI_Status status;
	//MPI_Request request;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int tag1 = 1, tag2 = 2, tag3 = 0;
	int isDone = 0, ii = 0, jj = 0;
	int Stolb_blokov_otprv = 0;
	int polnih_blokov_otpr = 0;
	int proshlo_str_blokov = 0;
	
	if (n < size)
	{
		limit = n + 1;// if kolichestvo processov bolwe razmera matrix
		N = 1;
		M = 1;
		Nn = 1;
		//is_n_menwe_size = true;
	}
	else
	{
		limit = size;  // if kolichestvo processov menwe ili ravno razmeru matrix
	// ceil - okruglenie do bolwogo
	// floor - okrug do menwego
		//N = ceil(n / limit-1); // kol. blokov strok

		//K = (n*m) - kol_stolb_v_posled_stolb_bloke    nihuya eto ne tak??// ili tAK?
		N = (n / (limit - 1.0)) + 1;
		M = (m / (limit - 1.0)) + 1;
		//Nn = floor((double)n / ((double)limit - 1));
		//M = ceil(m / limit-1); //kol. blok stolb
		//Mm = floor(m / limit - 1);
		Nn = N - 1;
		Mm = M - 1;
		kol_stolb_v_posled_stolb_bloke = m - Mm*(limit-1);
		kol_strok_v_posled_str_bloke = n - Nn*(limit-1);
		//is_n_menwe_size = false;

	}
	if (rank == 0) {
		//cout << kol_stolb_v_posled_stolb_bloke<<endl;
	   // cout << kol_strok_v_posled_str_bloke<<endl;
		cout << N << endl;
		cout << Nn << endl;
		FillMatrix(A, B);
		//FillVector(v);
		Zapis_v_File();
		//read_Vector();
		read_Matrix();

		for (int end_1 = 0; end_1 < (n * n); end_1 += end_1_otprav)  // 4et  poka viglyadit nepravilno
		{
			end_1_otprav = 0;
			//cout << "some part of all ranks" << " started with " << end_1 << endl;
			//cout << " poslednee ii - " << ii << endl;
			for (int i = 0; i < limit - 1; i++)
			{
				if ((((limit - 1) * proshlo_str_blokov) + i) < n)
				{
					for (int j = 0; j < limit - 1; j++)
					{
						//cout << i << "ne startuet" << endl;
						if ((((limit - 1) * Stolb_blokov_otprv) + j) < n)
						{
							//cout << j << "ne startuet" << endl;
							ii = (((limit - 1) * proshlo_str_blokov) + i); // Need to make actual schetchik stroki.
							jj = (((limit - 1) * Stolb_blokov_otprv) + j);
							//cout << "Otpr " << ii << " stroku, " << jj << " stolbec " << endl;
							vzat_vector_iz_matrix(A1, ii, B1, jj/*,1*/);

							MPI_Ssend(&k, m, MPI_DOUBLE, j + 1, tag1, MPI_COMM_WORLD);

							MPI_Rsend(&l, m, MPI_DOUBLE, j + 1, tag2, MPI_COMM_WORLD);
							//cout << "Otpr 1 rannk = " << rank << " " << ii << " stroku, " << jj << " stolbec " << endl;

							Mesto[0] = ii;
							Mesto[1] = jj;
							//cout << " otrpavlau " << ii << " stroku " << jj << " stolbec na " << j + 1 << " proccess" << endl;
							MPI_Rsend(&Mesto, 2, MPI_INT, j + 1, 5, MPI_COMM_WORLD);
							fflush(stdout);
							end_1_otprav++;
						}
						else
						{
							j = n * n;
						}
					}
				}
				else
				{
					i = (n * n);

				}
			}
	
			Stolb_blokov_otprv++;
			if (jj + 1 == n)
			{
				//ii += limit - 1;
				//jj = 0;
				//polnih_blokov_otpr++;
				Stolb_blokov_otprv = 0;
				proshlo_str_blokov++;
			}
		/*	else
			{
				ii = proshlo_str_blokov * (limit - 1);

			}*/
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				// tag3 = j;
				MPI_Recv(&(Temp[0]), 3, MPI_DOUBLE, MPI_ANY_SOURCE, tag3, MPI_COMM_WORLD, &status);
				if (Temp[1] < n && Temp[2] < n)
				{
					h = Temp[1];
					g = Temp[2];
					C[h][g] = Temp[0];
				}
				else
				{
					j--;
				}
				//cout << " gavno rabotai " << Temp[0] << endl;
				//cout << "Priem 2 real stroka = " << h <<" real stolbec "<< g << " C = " << C[h][g] << endl;
				fflush(stdout);
			}
		}

		Zapix_otvetov_v_File(C/*,d*/);
	}
	if ((rank > 0) && (rank <= kol_strok_v_posled_str_bloke)) {
		//for (int end = 0; end < m; end += end_1_priem)  // 4et  poka viglyadit nepravilno
		//{
		//	end_1_priem = 0;
		MPI_Request request;
		for (int j = 0; j < N*n; j++)
		{
			//tag3 = 0;
			//cout << "pppppppppPPPPPPPPPPAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA____________APPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPphui_PPPPPPPPPPPPPPPPPPPPPPP" << endl;
			fflush(stdout);
			MPI_Recv(&(k[0]), m, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);
			MPI_Recv(&(l[0]), m, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);
			MPI_Recv(&(Mesto[0]), 2, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);

			//cout << " prinal " << Mesto[0] << " stroku " << Mesto[1] << " stolbec na " << rank << " proccess" << endl;
			Temp[0] = peremnoj_vector_na_vector(k, l);
			Temp[1] = Mesto[0];
			Temp[2] = Mesto[1];
			//cout << " cam otvet = " << Temp[0] << "; ego stroka = " << Temp[1] << "; Ego stolbec = " <<Temp[2]<< endl;
			MPI_Isend(&Temp, 3, MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD, &request);
			MPI_Request_free(&request);

			//		end_1_priem++;
		}
		//	}

		//cout << rank << " all counting is done. Writing the answers" << endl;
	}
	
	//if (is_n_menwe_size)   // nahuya a glavnoe za4em? Ono zhe i beZ etogo rabotaet, da i s nim rabotaet, tak ono je ne nado, za4em ono togda rabotaet??
	//{
	//	limit_1 = limit;
	//}
	//else
	//{
	//	limit_1 = limit + 1;
	//}
	
	if ((rank > kol_strok_v_posled_str_bloke) && (rank < limit)) {
		double buff[1000];
		MPI_Buffer_attach(&buff, n* m * sizeof(double));    //rugaetsa na buff tipo ispolzovanie neinitialized memory kekl
		for (int j = 0; j < Nn*n; j++)
		{
			//tag3 = 0;
			MPI_Recv(&(k[0]), m, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);
			MPI_Recv(&(l[0]), m, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);
			MPI_Recv(&(Mesto[0]), 2, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);

			//cout << " prinal " << Mesto[0] << " stroku " << Mesto[1] << " stolbec na " << rank << " proccess" << endl;
			Temp[0] = peremnoj_vector_na_vector(k, l);
			Temp[1] = Mesto[0];
			Temp[2] = Mesto[1];
			//cout << " cam otvet = " << Temp[0] << "; ego stroka = " << Temp[1] << "; Ego stolbec = " << Temp[2] << endl;
			fflush(stdout);
			MPI_Bsend(&Temp, 3, MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD);
			//cout << rank << " second of all counting is done. Writing the answers" << endl;
			//		end_1_priem++;
		}
	}
	endtime = MPI_Wtime();
	printf("vipolnenie zanyalo %f seconds\n", endtime - starttime);
	MPI_Finalize();
	//MPI_Finalize();
	return 1;
};
