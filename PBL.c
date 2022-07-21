#include<stdio.h>
#include<windows.h>
#include <time.h>
typedef struct SquareMatrix{
	double e[100][100];
	int size;
}matrix;
typedef struct Vector{
	double e[100];
	int size;
}vector;
int getData(matrix* A, vector* B)//Đọc dữ liệu từ file
{
	char path[100];
	printf("Enter file name: ");
	scanf("%s",&path);//nhập tên file cần nhập
	FILE* f=fopen(path,"r");//mở file để đọc
	if (f!=NULL){//nếu tìm thấy file sẽ tiếp tục ngược lại trả về 0
		char c;
		A->size=0;   B->size=0; // khỏi tạo kích thước ma trận và vector về 0
		fscanf(f,"%d",&(A->size));// đọc kí tự đầu tiên là kích thước ma trận và veotor
		B->size=A->size;
		if(A->size!=0){//nếu kích thước là 0 thì file nhập vào là rỗng
			for (int i=0;i<A->size;i++)
			{

				for (int j=0;j<A->size;j++) //nhập n kí tự đầu tiên dòng i cho ma trận A[i][j]
					if(fscanf(f,"%lf",&(A->e[i][j]))==EOF) { //nếu như không tìm thấy kí tự nhập vào thì dữ liệu sai
						printf("Invalid %s file data!\n(Please enter again file data)\n",path);
						fclose(f);
						return 0;
					}
				if(fscanf(f,"%lf",&(B->e[i]))==EOF) { //nhập phần tử n+1 cho vector B[i], nếu như không tìm thấy kí tự nhập vào thì dữ liệu sai
						printf("Invalid %s file data!\n(Please enter again file data)\n",path);
						fclose(f);
						return 0;
					}
			}
			printf("Read file %s successfully!\n",path);
		}
		else{//file nhập vào là file rỗng lỗi đọc file trả về 0
			printf("The file %s is empty!\n",path);
			printf("**********************************************************\n");
			fclose(f);return 0;
		}
	}
	else {//Không tìm thấy file trả về 0
		printf("File %s not found!\n",path);
		printf("**********************************************************\n");
		fclose(f);return 0;
	}
	printf("**********************************************************\n");
	fclose(f);

	return 1;//Đọc file thành công trả về 1
}

void writeSolution(matrix A,vector X,vector B){
	char path[100];
	printf("Enter file name: ");
	scanf("%s",&path);
	FILE* f=fopen(path,"a");
	fprintf(f,"\n");
	printf("Write file %s successfully!\n",path);
	for (int i=0;i<A.size;i++){
		for (int j=0;j<A.size;j++)
			fprintf(f,"%8.3lf ",A.e[i][j]);
		fprintf(f,"%8.3lf ",B.e[i]);
		fprintf(f,"\n");
	}
	fprintf(f,"Solution of the system of linear equations X=( ");
	for (int i=0;i<X.size-1;i++) fprintf(f,"%.3lf, ",X.e[i]);
	fprintf(f,"%.3lf)\n",X.e[X.size-1]);
	printf("**********************************************************\n");
	fclose(f);
}
void displayProblem(matrix A,vector B)//in ra ma tran
{
	for (int i=0;i<A.size;i++)
	{
		for (int j=0;j<A.size;j++)
			printf("  %8.3lf",A.e[i][j]);
		printf("  %8.3lf",B.e[i]);
		printf("\n");
	}
}
void displayMatrix(matrix A)//in ra ma tran
{
	for (int i=0;i<A.size;i++)
	{
		for (int j=0;j<A.size;j++)
			printf("  %8.3lf",A.e[i][j]);
		printf("\n");
	}
}
void displayResult(vector X){
	printf("Solution of the system of linear equations X=( ");
	for (int i=0;i<X.size-1;i++) printf("%.3lf, ",X.e[i]);
	printf("%.3lf)\n",X.e[X.size-1]);
}
void swap(double* x, double* y)// hoan vi 2 gia tri x va y
{
	double temp=*x;
	*x=*y;
	*y=temp;
}
double det(matrix A)//tinh dinh thuc cua ma tran A la detA (|A|)
{
	double detA=1,k;
	for (int i=0;i<A.size-1;i++)
	{
		for (int j=i+1;j<A.size;j++)
		if (A.e[j][i]!=0 && A.e[i][i]!=0)
		{
			k=A.e[j][i]/A.e[i][i];
			for (int col=i;col<A.size;col++)
				A.e[j][col]-=A.e[i][col]*k;
		}
		else if (A.e[i][i]==0)//khong the dung gauss phai hoan vi
			for (int k=i+1;k<A.size;k++)
				if (A.e[k][i]!=0)//hoan vi 2 hang cua ma tran
				{
					for (int col=0;col<A.size;col++)
						swap(&A.e[k][col],&A.e[i][col]);
					j--;
					detA*=-1;
					break;
				}
		detA*=A.e[i][i];
	}
	return detA*A.e[A.size-1][A.size-1];
}
vector solveByGauss(matrix A,vector B){
	double k;   vector X;
	X.size=B.size;
	for (int i=0;i<A.size-1;i++)
	{
		for (int j=i+1;j<A.size;j++)
		if (A.e[j][i]!=0 && A.e[i][i]!=0)
		{
			k=A.e[j][i]/A.e[i][i];
			for (int col=i;col<A.size;col++)
				A.e[j][col]-=A.e[i][col]*k;
			B.e[j]-=B.e[i]*k;
			if(k>0) printf("Transformation:    R%d - %.3lf x R%d----> R%d\n",j+1,k,i+1,j+1);
			else printf("Transformation:    R%d + %.3lf x R%d----> R%d\n",j+1,k*-1,i+1,j+1);
			displayProblem(A,B);
			printf("\n");
		}
		else if (A.e[i][i]==0){//khong the dung gauss phai hoan vi
			for (int k=i+1;k<A.size;k++)
				if (A.e[k][i]!=0)//hoan vi 2 hang cua ma tran
				{
					for (int col=0;col<A.size;col++)
						swap(&A.e[k][col],&A.e[i][col]);
					swap(&B.e[k],&B.e[i]);
					j--;
					printf("Transformation:    R%d <----> R%d\n",i+1,k+1);
					displayProblem(A,B);
					printf("\n");
					break;
				}
		}
		if(!A.e[i][i]){
			X.size=0;
			return X;
		}
	}
	if(!A.e[A.size-1][A.size-1]){
			X.size=0;
			return X;
		}
	for(int i=A.size-1;i>=0;i--){
		X.e[i]=B.e[i];
		for (int j=A.size-1;j>i;j--) X.e[i]-=A.e[i][j]*X.e[j];
		X.e[i]/=A.e[i][i];
	}
	return X;
}
vector Solution(matrix A,vector B){//tương tự gauss dùng để tính thời gian
	double k;   vector X;
	X.size=B.size;
	for (int i=0;i<A.size-1;i++)
	{
		for (int j=i+1;j<A.size;j++)
		if (A.e[j][i]!=0 && A.e[i][i]!=0)
		{
			k=A.e[j][i]/A.e[i][i];
			for (int col=i;col<A.size;col++)
				A.e[j][col]-=A.e[i][col]*k;
			B.e[j]-=B.e[i]*k;
		}
		else if (A.e[i][i]==0){//khong the dung gauss phai hoan vi
			for (int k=i+1;k<A.size;k++)
				if (A.e[k][i]!=0)//hoan vi 2 hang cua ma tran
				{
					for (int col=0;col<A.size;col++)
						swap(&A.e[k][col],&A.e[i][col]);
					swap(&B.e[k],&B.e[i]);
					j--;
					break;
				}
		}
		if(!A.e[i][i]){
			X.size=0;
			return X;
		}
	}
	if(!A.e[A.size-1][A.size-1]){
			X.size=0;
			return X;
		}
	for(int i=A.size-1;i>=0;i--){
		X.e[i]=B.e[i];
		for (int j=A.size-1;j>i;j--) X.e[i]-=A.e[i][j]*X.e[j];
		X.e[i]/=A.e[i][i];
	}
	return X;
}
vector Solution2(matrix A,vector B){//tương tự cramer dùng để tính thời gian
	matrix M;   	 vector X;
	M.size=A.size;
	X.size=B.size;
	double detM,detA=det(A);
	if(detA!=0){
		for(int k=0;k<A.size; k++){
			for (int i=0;i<A.size;i++)
			{
				for (int j=0;j<k;j++) M.e[i][j]=A.e[i][j];

				for (int j=k+1;j<A.size;j++) M.e[i][j]=A.e[i][j];

				M.e[i][k]=B.e[i];
			}
			detM=det(M);
			X.e[k]=detM/detA;
		}
	}
	else{
		X.size=0;
	}
	return X;
}
vector solveByCrame(matrix A,vector B){
	matrix M;   	 vector X;
	M.size=A.size;
	X.size=B.size;
	double detM,detA=det(A);
	printf("Matrix A:\n");
	displayMatrix(A);
	printf("\ndet(A)=%.3Lf\n\n",detA);
	if(detA!=0){
		for(int k=0;k<A.size; k++){
			for (int i=0;i<A.size;i++)
			{
				for (int j=0;j<k;j++) M.e[i][j]=A.e[i][j];

				for (int j=k+1;j<A.size;j++) M.e[i][j]=A.e[i][j];

				M.e[i][k]=B.e[i];
			}
			detM=det(M);
			X.e[k]=detM/detA;
			printf("Matrix A%d:\n",k+1);
			displayMatrix(M);
			printf("\ndet(A%d)=%.3Lf---->X%d=det(A%d)/det(A)=%0.3Lf\n\n",k+1,detM,k+1,k+1,X.e[k]);
		}
	}
	else{
		X.size=0;
	}
	return X;
}

int main()
{
	clock_t start, end;
	int option,Data=0;
	char ok;
	matrix A;
	vector B,X;
	do
	{
	printf("---------------------Function table---------------------\n");
	printf("\t1.Read matrix from input file\n");
	printf("\t2.Display matrix from file\n");
	printf("\t3.Show solving steps by Crame\n");
	printf("\t4.Show solving steps by Gauss\n");
	printf("\t5.Display the solution\n");
	printf("\t6.Write the solution of the equation in the file\n");
	printf("\t7.Compare the execution time of gauss and crame\n");
	printf("\t8.Exit\n");
	printf("\n****Please enter options: ");
	scanf("%d",&option);
	printf("\n");
	switch (option)
		{
			case 1:
				Data=getData(&A,&B);
				break;
			case 2:
				if(!Data)printf("Data is not available\n"),Data=getData(&A,&B);
				if(!Data) break;
				displayProblem(A,B);
				break;
			case 3:
				if(!Data)printf("Data is not available\n"),Data=getData(&A,&B);
				if(!Data) break;
					X=solveByCrame(A,B);
					if(X.size) displayResult(X);
					else printf("No Solution!\n");
				break;
			case 4:
				if(!Data)printf("Data is not available\n"),Data=getData(&A,&B);
				if(!Data) break;
					X=solveByGauss(A,B);
					if(X.size) displayResult(X);
					else printf("No Solution!\n");
				break;
			case 5:
				if(!Data)printf("Data is not available\n"),Data=getData(&A,&B);
				if(!Data) break;
					X=Solution(A,B);
					if(X.size) displayResult(X);
					else printf("No Solution!\n");
				break;
			case 6:
				if(!Data)printf("Data is not available\n"),Data=getData(&A,&B);
				if(!Data) break;
					X=Solution(A,B);
					if(X.size) writeSolution(A,X,B);
					else printf("No Solution!\n");
				break;
			case 7:
				if(!Data)printf("Data is not available\n"),Data=getData(&A,&B);

				start=clock();
				X=Solution(A,B);
				end	=clock();
				printf("Ecution time of Gauss: %lf\n",(double)(end-start)/CLOCKS_PER_SEC);
				start=clock();
				X=Solution2(A,B);
				end	=clock();
				printf("Ecution time of Crame: %lf\n",(double)(end-start)/CLOCKS_PER_SEC);
				break;
			case 8:
				break;
			default:printf("This function could not be found!\n");
		}
		if(option==8) break;
		else
		printf("\nDo you want to continue?");
		do{
			printf("Yes or No enter (y/n): ");
			fflush(stdin);
			scanf("%c",&ok);
			if(ok=='y' || ok=='Y') system("cls");
			else if(ok=='n' || ok=='N')option=0;
		}while(ok!='y'&&ok!='Y'&&ok!='n'&&ok!='N');

	}while (option);
	return 0;
}

