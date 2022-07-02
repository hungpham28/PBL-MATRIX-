#include<stdio.h>
void enter_matrix(int r,int c,double A[100][100])//nhap ma tran
{
	for (int i=0;i<r;i++)
		for (int j=0;j<c;j++)
			scanf("%lf",&A[i][j]);
}
void print_matrix(int r,int c, double A[100][100])//in ra ma tran
{
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
			printf("  %0.4lf",A[i][j]);
		printf("\n"); 
	}
}
void swap(double* x, double* y)// hoan vi 2 gia tri x va y
{
	double temp=*x;
	*x=*y;
	*y=temp;
}
void multi_matrix(int r1,int c1, double A[100][100],int r2, int c2, double B[100][100], double C[100][100] )//phep nhan ma tran AxB
{
	for (int i=0;i<r1; i++)
	for (int j=0;j<c2;j++)
	{
	 	C[i][j]=0;
	 	for (int rc=0;rc<c1;rc++)
	 		C[i][j]+=A[i][rc]*B[rc][j];
	}
}
void add_matrix( int r,int c, double A[100][100],double B[100][100], double C[100][100])//phep cong ma tran A+B
{
	for (int i=0;i<r; i++)
	for (int j=0;j<c;j++)
		C[i][j]=A[i][j]+B[i][j];
}
double det(int n, double A[100][100])//tinh dinh thuc cua ma tran A la detA (|A|)
{
/*--------------------Phuong phap Gauss-------------------*/
	double detA=1,M[100][100],k;
	for (int i=0;i<n; i++)//tao 1 ma tran giong A
		for (int j=0;j<n;j++)
			M[i][j]=A[i][j];
//			print_matrix(n,n,M);
//			printf("\n");
	for (int i=0;i<n-1;i++)
	{
		for (int j=i+1;j<n;j++)
		if (M[j][i]!=0 && M[i][i]!=0)
		{
			k=M[j][i]/M[i][i];
			for (int col=i;col<n;col++)
				M[j][col]-=M[i][col]*k;
//			print_matrix(n,n,M);
//			printf("\n");
		}
		else if (M[i][i]==0)//khong the dung gauss phai hoan vi
			for (int k=i+1;k<n;k++)
				if (M[k][i]!=0)//hoan vi 2 hang cua ma tran
				{
					for (int col=0;col<n;col++) 
						swap(&M[k][col],&M[i][col]);
					j--;
					detA*=-1;
//			print_matrix(n,n,M);
//			printf("------------------------------------\n");
					break;
				}
		detA*=M[i][i];
	}
	return detA*M[n-1][n-1];
}
void transpose(int r,int c, double A[100][100],double tranA[100][100])//ma tran chuyen vi cua A
{
	for (int i=0;i<c; i++)
		for (int j=0;j<r;j++)
			tranA[i][j]=A[j][i];
}

void delete_col(int r,int c, double A[100][100], int col)//xoa cot (col) ma tran A
{
	for (int i=col;i<c;i++)
		for (int j=0;j<r;j++)
			A[j][i-1]=A[j][i];
}
void delete_row(int r,int c, double A[100][100], int row)//xoa hang (row) ma tran A
{
	for (int i=row;i<r;i++)
		for (int j=0;j<c;j++)
			A[i-1][j]=A[i][j];
}

int inverse_matrix(int n,double A[100][100],double I[100][100])//ma tran nghich dao
{
/*--------------------Phuong phap Gauss-------------------*/

	double M[100][100],k;
	for (int i=0;i<n; i++)
	{	//khoi tao ma tran don vi I
		for (int j=0;j<n;j++)
			I[i][j]=0;//I la ma tran don vi 1
		I[i][i]=1;
	}

	for (int i=0;i<n; i++)//tao 1 ma tran giong A
		for (int j=0;j<n;j++)
			M[i][j]=A[i][j];
			
	for (int i=0;i<n-1;i++)
	{
		for (int j=i+1;j<n;j++)
		if (M[j][i]!=0 && M[i][i]!=0)
		{
			k=M[j][i]/M[i][i];
			for (int col=0;col<n;col++)
				M[j][col]-=M[i][col]*k , I[j][col]-=I[i][col]*k;
//			print_matrix(n,n,M);
//			printf("matrix A\n\n");
//			print_matrix(n,n,I);
//			printf("matrix I\n------------------------------------\n");
		}
		else if (M[i][i]==0)//khong the dung gauss phai hoan vi
			for (int k=i+1;k<n;k++)
				if (M[k][i]!=0)//hoan vi 2 hang cua ma tran
				{
					for (int col=0;col<n;col++) 
						swap(&M[k][col],&M[i][col]),swap(&I[k][col],&I[i][col]);
					j--;
//					print_matrix(n,n,M);
//					printf("matrix A\n\n");
//					print_matrix(n,n,I);
//					printf("matrix I\n------------------------------------\n");
					break;
				}
	}
	for (int i=0;i<n;i++) if (M[i][i]==0) return 0;
	for (int i=0;i<n;i++)// dua duong cheo chinh ve 1
	{
		k=M[i][i];
		for (int j=0;j<n;j++)
			M[i][j]=M[i][j]/k,I[i][j]=I[i][j]/k;
//		print_matrix(n,n,M);
//		printf("matrix A\n\n");
//		print_matrix(n,n,I);
//		printf("matrix I\n------------------------------------\n");
		
	}
	for (int i=n-1;i>0;i--)
	{
		for (int j=i-1;j>=0;j--)
		if (M[j][i]!=0 && M[i][i]!=0)
		{
			k=M[j][i];
			for (int col=n-1;col>=0;col--)
				M[j][col]-=M[i][col]*k , I[j][col]-=I[i][col]*k;
//			print_matrix(n,n,M);
//			printf("matrix A\n\n");
//			print_matrix(n,n,I);
//			printf("matrix I\n------------------------------------\n");
		}
	}
	return 1;
}

int main()
{
	int n,option;
	printf("-------Bang chuc nang cac phep toan voi ma tran-------\n");
	printf("\t0.Nhan phim 0 de thoat\n\t1.Phep cong ma tran A+B\n\t2.Xoa cot ma tran\n\t3.Xoa hang ma tran\n\t4.Phep nhan ma tran AxB\n\t5.Ma tran chuyen vi cua A\n\t6.Tinh dinh thuc det(A)\n\t7.Ma tran nghich dao cua A");
	do
	{
	printf("\n****Moi ban chon phep tinh voi ma tran: ");
	scanf("%d",&option);
	switch (option)
		{
			case 0:break;
			case 1:{
				int m,n;
				printf("Nhap kich thuoc mxn cho 2 ma tran A va B: ");
				scanf("%d%d",&m,&n);
				double A[100][100],B[100][100],C[100][100];
				printf("Nhap ma tran A:\n");
				enter_matrix(m,n,A);
				printf("Nhap ma tran B:\n");
				enter_matrix(m,n,B);
				printf("Phep cong ma tran A+B la\n");
				add_matrix(m,n,A,B,C);
				print_matrix(m,n,C);
				break;
			}
			case 2:{
				int m,n,col;
				printf("Nhap kich thuoc mxn cho ma tran A: ");
				scanf("%d%d",&m,&n);
				double A[100][100];
				printf("Nhap ma tran A:\n");
				enter_matrix(m,n,A);
				printf("Nhap cot muon xoa:");
				scanf("%d",&col);
				delete_col(m,n,A,col);
				print_matrix(m,--n,A);
				break;
			}
			case 3:{
				int m,n,row;
				printf("Nhap kich thuoc mxn cho ma tran A: ");
				scanf("%d%d",&m,&n);
				double A[100][100];
				printf("Nhap ma tran A:\n");
				enter_matrix(m,n,A);
				printf("Nhap hang muon xoa:");
				scanf("%d",&row);
				delete_row(m,n,A,row);
				print_matrix(--m,n,A);
				break;
			}
			case 4:{
				int m1,n1,m2,n2;
				do{
				printf("Nhap kich thuoc mxn cho ma tran A: ");
				scanf("%d%d",&m1,&n1);
				printf("Nhap kich thuoc mxn cho ma tran B: ");
				scanf("%d%d",&m2,&n2);
				if (n1!=m2) printf("so cot ma tran A phai = so hang ma tran B!\n");
				}while(n1!=m2);
				double A[100][100],B[100][100],C[100][100];
				printf("Nhap ma tran A:\n");
				enter_matrix(m1,n1,A);
				printf("Nhap ma tran B:\n");
				enter_matrix(m2,n2,B);
				printf("Phep nhan ma tran AxB la\n");
				multi_matrix(m1,n1,A,m2,n2,B,C);
				print_matrix(m1,n2,C);
				break;
			}
			case 5:{
				int m,n;
				printf("Nhap kich thuoc mxn cho ma tran A : ");
				scanf("%d%d",&m,&n);
				double A[100][100],tranA[100][100];
				printf("Nhap ma tran A:\n");
				enter_matrix(m,n,A);
				printf("Ma tran chuyen vi cua A la\n");
				transpose(m,n,A,tranA);
				print_matrix(n,m,tranA);
				break;
			}
			case 6:{
				int m,n;
				printf("Nhap kich thuoc n cho ma tran vuong A : ");
				scanf("%d",&n);
				double A[100][100];
				printf("Nhap ma tran A:\n");
				enter_matrix(n,n,A);
				printf("Dinh thuc cua ma tran det(A)=%0.4lf ",det(n,A));
				break;
			}
			case 7:{
				int m,n;
				printf("Nhap kich thuoc n cho ma tran vuong A : ");
				scanf("%d",&n);
				double A[100][100],inverA[100][100];
				printf("Nhap ma tran A:\n");
				enter_matrix(n,n,A);
				if(inverse_matrix(n,A,inverA)!=0)
				printf("Ma tran nghich dao cua A: \n"),
				print_matrix(n,n,inverA);
				else printf("Dinh thuc Det(A)=0.Ma tran khong kha nghich!")	;			
				break;
			}
			default:printf("Khong co chuc nang nay!");
		}
	}while (option);

}
//---------------------EXAMPLE--------------------
/*Matrix 1
7 7 9
1 1 5
7 6 8
Matrix 2
1 2 8 7 6
7 8 6 7 5
4 5 7 8 6
4 8 6 2 1

-4 -4 0 2 2 
-5 -2 2 3 5 
2 -2 4 1 1
4 4 -2 4 -2
Matrix 3
1 2 4 2
2 4 2 1
1 2 4 7
matrix 4
1 2
4 5
1 2
7 8

1 -4 
2 1
4 5 
5 -5
matrix 5
0 -3 -2 3 -1 0 0 -4
4 0 -2 -4 -2 -5 1 -5
-2 -5 -2 -3 -2 1 -3 -3
-4 -5 -4 -1 3 -3 2 -2
4 -4 -1 2 -5 5 4 -1
-1 -1 0 -2 -2 -3 5 3
4 5 1 1 2 -4 4 -4
-2 0 0 1 2 -5 -2 5
matrix 6
0 -3 -5 -3 1 -3 5 4 -3 0 
-5 3 0 2 1 -5 -2 0 -5 4 
-2 -4 -5 3 -3 -2 5 -3 2 3
1 -2 0 -2 1 0 5 -5 3 4
-5 -4 1 4 4 1 2 3 4 -3
-3 3 -4 2 0 4 0 -4 -3 3
-2 4 4 -3 -3 -4 -1 -2 -2 -4
3 -4 1 -2 -1 5 5 3 0 1
5 1 -3 -5 -5 0 -5 -1 -4 -1
-1 -2 5 -1 2 -3 -5 0 3 0
matrix 7
-5 0 5 -1 -5 -3 -5 -2 -1 3 3 -1 4 0 2 -3 -1 5 -4 -5 
1 2 2 4 1 0 1 -4 -1 -5 1 2 3 5 -2 -3 4 4 -5 -5 
5 -3 1 -1 -2 2 2 1 3 -4 5 -4 2 2 2 -4 -3 2 -5 -3
0 5 1 4 -5 5 -3 2 2 1 0 5 1 3 5 1 1 -4 1 4
-4 -5 -4 4 3 0 -4 -2 -4 1 -4 0 3 3 -4 -4 -3 -1 -3 2
3 -1 3 3 -1 2 4 -1 5 0 -3 -4 -3 -5 -2 1 -1 -1 -4 1
2 4 4 0 -3 0 2 -2 -5 5 -4 4 4 -1 4 0 -5 2 2 -5
-5 -2 -3 -4 -4 -4 2 -5 -4 1 3 -5 -5 5 -2 0 -3 -1 -5 -1
-2 -1 -2 -2 5 1 3 -4 5 5 0 -1 -3 -3 3 2 3 3 3 1
1 -3 -2 -5 3 -2 -1 -1 -5 5 -2 0 0 -1 5 -1 3 2 3 5
4 -2 2 2 2 5 -3 2 -5 -1 -4 -4 0 -3 2 -1 -4 1 3 0
0 4 -3 -4 -3 4 5 -2 -2 1 0 3 1 2 -1 0 -4 1 -4 0
-4 5 5 0 3 -5 5 5 -5 -3 3 -4 -1 5 5 5 2 0 -5 -4
-1 3 4 1 -1 1 0 0 0 -1 3 2 2 0 0 -5 5 2 4 0
1 5 1 -2 -5 1 -2 3 -1 -4 -3 1 -5 -3 4 0 0 -3 -5 0
2 -4 0 1 -3 -4 -3 1 4 -1 0 4 1 -4 2 2 5 -5 -1 2
4 4 -5 -4 4 3 1 5 0 2 -2 2 -3 5 2 -3 0 4 4 0
4 3 5 -4 0 -5 2 0 3 4 2 -4 0 2 3 -1 -2 -5 -5 -2
-1 -3 1 5 0 4 5 0 -2 -4 -2 -3 -2 2 4 -5 -1 3 -2 3
-1 -2 1 -1 0 -2 4 0 5 5 -4 -1 -2 0 1 -4 4 -5 3 0
*/
