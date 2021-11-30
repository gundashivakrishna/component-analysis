#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define MROW 800

int N;
int no_rows[3];
bool diff_datasets;
///////////////////////////////////////


void getCofactor(float A[][N], float temp[][N], int p, int q, int n) 
{ 
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q) 
            { 
                temp[i][j++] = A[row][col]; 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 

int determinant(float A[][N], int n) 
{ 
    int D = 0; // Initialize result 
  
    //  Base case : if matrix contains single element 
    if (n == 1) 
        return A[0][0]; 
  
    float temp[N][N]; // To store cofactors 
  
    int sign = 1;  // To store sign multiplier 
  
     // Iterate for each element of first row 
    for (int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        getCofactor(A, temp, 0, f, n); 
        D += sign * A[0][f] * determinant(temp, n - 1); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
  
    return D; 
} 
  
// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(float A[][N],float adj[][N]) 
{ 
    if (N == 1) 
    { 
        adj[0][0] = 1; 
        return; 
    } 
  
    // temp is used for storing cofactors of A[][] 
	int sign = 1;
	float temp[N][N]; 

	for (int i=0; i<N; i++) 
	{ 
		for (int j=0; j<N; j++) 
		{ 
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N);
            sign = ((i+j)%2==0)? 1: -1;
            adj[j][i] = (sign)*(determinant(temp, N-1)); 
        } 
    } 
} 
  

bool inverse(float A[][N], float inverse[][N]) 
{ 
    // Find determinant of A[][] 
    float det = determinant(A, N); 
    if (det == 0) 
    { 
        //printf("Singular matrix, can't find its inverse\n"); 
        return false; 
    } 
  
    // Find adjoint 
    float adj[N][N]; 
    adjoint(A, adj); 
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i=0; i<N; i++) 
        for (int j=0; j<N; j++) 
            inverse[i][j] = adj[i][j]/det; 
  
    return true; 
} 
//***********************



void transpose(float array[][N],int rows,int columns,float t[][N]){
	 for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j) {
            t[j][i] = array[i][j];
        }
	return ;
}

void multiplyMatrices(float first[][10],float second[][10],float mult[][10], int r1, int c1, int r2, int c2) {

    // Initializing elements of matrix mult to 0.
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            mult[i][j] = 0;
        }
    }

    // Multiplying first and second matrices and storing in mult.
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            for (int k = 0; k < c1; ++k) {
                mult[i][j] += first[i][k] * second[k][j];
            }
        }
    }
	return;
}


void calculate_C(float array[][N],float C[][N],int new_rows){
	float xi_minus_xj[N][N],xi_minus_xj_transpose[N][N];
	for(int i=0;i<new_rows;i++){
		//only unique combinations of i and j or not?
		for(int j=0;j<new_rows;j++){
			for(int k=0;k<N;k++){
				xi_minus_xj[k][k]=array[i][k]-array[j][k];
				transpose(xi_minus_xj,new_rows,N,xi_minus_xj_transpose);
				multiplyMatrices(xi_minus_xj,xi_minus_xj_transpose,C,N,N,N,N);
			}
		}
	}
	return;
}

//*************************


///////////////////////////////////////

bool isPresent(float class,float* classes,int sz){
	sz--;
	for (int i = 0; i < sz; ++i)
	{
		if(classes[i]==class) return true;
	}
	return false;
}


void eigen(float a[N][N],float* x)
{
	 float x_new[N];
	 float temp, lambda_new, lambda_old, error=0.00001;
	 int i,j,n=N, step=1;
	 lambda_old = 1;
	 /* Multiplication */
	 up:
	 for(i=0;i<n;i++)
	 {
		  temp = 0.0;
		  for(j=0;j<n;j++)
		  {
		   	temp = temp + a[i][j]*x[j];
		  }
		  x_new[i] = temp;
	 }
	 /* Replacing */
	 for(i=0;i<n;i++)
	 {
	  	x[i] = x_new[i];
	 }
	 /* Finding Largest */
	 lambda_new = fabs(x[0]);
	 for(i=1;i<n;i++)
	 {
		  if(fabs(x[i])>lambda_new)
		  {
		   	lambda_new = fabs(x[i]);
		  }
	 }
	 /* Normalization */
	 for(i=0;i<n;i++)
	 {
	  	x[i] = x[i]/lambda_new;
	 }

	 if(fabs(lambda_new-lambda_old)>error)
	 {
		  lambda_old=lambda_new;
		  step++;
		  goto up;
	 }

	/*temp = 0;
    	for (int i = 0; i < N; i++)
    	{
        	temp += (x[i]) * (x[i]);
    	}
    	temp = sqrt(temp);
    	for (int i = 0; i < N; i++)
    	{
        	x[i] = x[i] / temp;
    	}*/
	 
}

void pca(int no_classes,int no_fields, char* fname, char* fname1){
	N=no_fields-1;
	FILE *fp1, *fp2;
    float a[MROW][no_fields], c[N][N];
    float eigV[N], temp, mean[N];

    //taking dataset and storing in an array, normalizing it
    fp1 = fopen(fname, "r");
    fp2 = fopen(fname1, "w+");
    //printf("1111111111\n");
    for (int i = 0; i < N; i++)
    {
        eigV[i] = 1;
    }
    //printf("22222222222\n");
    /*while (!feof(fp1))
    {
        for (int i = 0; i < COUNT; i++)
        {
            for (int j = 0; j < no_fields; j++)
            {
                fscanf(fp1, "%f ", &a[i][j]);
            }
        }
    }*/
    int rows=0;
	while(!feof(fp1)){
		// if(a){

		// }
		for(int i=0;i<no_fields;i++){
			fscanf(fp1,"%f",&(a[rows][i]));
			
		}
		//printf("%d\n",rows);
		//printf("%f\n",a[rows][no_fields-1]);
		rows++;
	}
	//if(a[rows][no_fields-1]-floor(a[rows][no_fields-1])!=0) rows--;
	if(a[rows-2][no_fields-1] != a[rows-1][no_fields-1]) rows--;
	no_rows[0]=rows;
	//printf("data reading done!!\n");
	float b[rows][N];
	float res[rows];
    for (int i = 0; i < N; i++)
    {
        temp = 0;
        for (int j = 0; j < rows; j++)
        {
            temp += a[j][i];
        }
        mean[i] = temp / rows;
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < no_fields; j++)
        {
            b[i][j] = a[i][j] - mean[j];
        }
    }

    //Covariance Matrix
    for (int j1 = 0; j1 < N; j1++)
    {
        for (int j2 = j1; j2 < N; j2++)
        {
            temp = 0;
            for (int i = 0; i < rows; i++)
            {
                temp += b[i][j1] * b[i][j2];
            }
            c[j1][j2] = temp / rows;
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = N; j > i; j--)
        {
            c[j][i] = c[i][j];
        }
    }

    //finding eigV and normalizing it
    //printf("a");
    eigen(c, eigV);
    //printf("b");
    //Projecting dataset onto eigV
    for (int i = 0; i < rows; i++)
    {
        temp = 0.0;
        for (int j = 0; j < N; j++)
        {
            temp += b[i][j] * eigV[j];
        }
        res[i] = temp;
    }
    //printf("rows are %d\n",rows);
    //printing result to a file
    for (int i = 0; i < rows; i++)
    {
        fprintf(fp2, "%f ", res[i]);
        fprintf(fp2, "%f\n", a[i][no_fields-1]);
    }

    fclose(fp1);
    fclose(fp2);
}

void fld(int no_classes,int no_fields, char* fname, char* fname1){
	N=no_fields-1;
	FILE* fp=fopen(fname,"r");
	FILE* fp2=fopen(fname1,"w+");
	// printf("Is serial number present in data: \n");
	// int a;
	// scanf("%d",&a);
	int tmp=0;
	float data[MROW][no_fields];
	float classes[no_classes];
	int nos[no_classes];
	memset(nos,0,sizeof(nos));
	int rows=0;
	while(!feof(fp)){
		// if(a){

		// }
		for(int i=0;i<no_fields;i++){
			fscanf(fp,"%f",&data[rows][i]);
	
		}
		//printf("%d\n",rows);
		//printf("%f\n",data[rows][no_fields-1]);
		rows++;
	}
	if((data[rows-1][0]||data[rows-1][1])==0){ rows-=1;}
	printf("Data Set reading done\n");
	//printf("No. of data points: %d\n",rows);
	no_rows[1]=rows;

	for(int i=0;i<rows;i++){
		if(!isPresent(data[i][no_classes-1],classes,tmp)){
			classes[tmp]=data[i][no_classes-1];
			tmp++;
		}
		if(tmp==no_classes) break;
	}
	tmp=0;

	float sum[no_classes][no_fields-1];
	memset(sum,0,sizeof(sum));
	for (int i = 0; i < no_classes; ++i)
	{
		for (int j = 0; j<rows; ++j)
		{
			if(data[j][no_classes-1]==classes[i]){
				nos[i]++;
				for (int k = 0; k < no_fields-1; ++k)
				{
					sum[i][k]+=data[j][k];
				}
			}
		}
	}
	float mean[no_classes][no_fields-1];
	for (int i = 0; i < no_classes; ++i)
	{
		for (int j = 0; j < no_fields-1; ++j)
		{
			mean[i][j]=(float)sum[i][j]/nos[i];
		}
	}

	//printf("%f\n",mean[0][0]);

	float Sw[no_fields-1][no_fields-1];
	memset(Sw,0,sizeof(Sw));
	for (int i = 0; i < no_classes; ++i)
	{
		for (int j = 0; j<rows; ++j)
		{
			if(data[j][no_classes-1]==classes[i]){
				for (int k = 0; k < no_fields-1; ++k)
				{
					for (int l = 0; l < no_fields-1; ++l)
					{
						Sw[k][l]+=(data[j][k]-mean[i][k])*(data[j][l]-mean[i][l]);
					}
				}

			}
		}
	}
	printf("Sw Computation done!!\n");
	float Sb[no_fields-1][no_fields-1];
	memset(Sb,0,sizeof(Sb));
	float mean_of_means[no_fields-1];
	memset(mean_of_means,0,sizeof(mean_of_means));

	for (int i = 0; i < no_fields-1; ++i){
		for (int j = 0; j < no_classes; ++j){
			mean_of_means[i]+=mean[j][i];
		}
		mean_of_means[i]=mean_of_means[i]/no_classes;
	}

	for(int i=0;i<no_classes;i++){
		for (int k = 0; k < no_fields-1; ++k)
		{
			for (int l = 0; l < no_fields-1; ++l)
			{
				Sb[k][l]+=nos[i]*((mean[i][k]-mean_of_means[k])*(mean[i][l]-mean_of_means[l]));
			}
		}
	}
	printf("Sb Computation done!!\n");
	float Sw_inv[no_fields-1][no_fields-1];
	inverse(Sw,Sw_inv);
	printf("Sw_inv Computation done!!\n");
	float Prod[no_fields-1][no_fields-1];
	memset(Prod,0,sizeof(Prod));
	for(int i=0;i<no_fields-1;i++){
		for(int j=0;j<no_fields-1;j++){
			for(int k=0;k<no_fields-1;k++){
				Prod[i][j]+=(Sw_inv[i][k]*Sb[k][j]);
			}
		}
	}

	float eigenvec[no_fields-1];
	memset(eigenvec,1,sizeof(eigenvec));
	eigen(Prod,eigenvec);
	printf("eigenvec Cmoputation done!!\n");
	float final_data[MROW][2];
	memset(final_data,0,sizeof(final_data));
	for(int i=0;i<rows;i++){
		for(int j=0;j<no_fields-1;j++){
			final_data[i][0]+=eigenvec[j]*data[i][j];
		}
		final_data[i][1]=data[i][no_fields-1];
	}
	//FILE* fdat=fopen(fname1,"w+");
	for(int i=0;i<rows;i++){
		fprintf(fp2,"%f  %f",(final_data[i][0]*1000),final_data[i][1]);
		if(i!=rows-1) fprintf(fp2,"\n");
	}
	printf("Final Data written to file %s\n",fname1);
	fclose(fp);
	fclose(fp2);
}



void metric_learning(int no_classes, int no_fields, char* fname, char* fname1){
   

	float sum_S[20][20],sum_D[20][20],sigma_S[20][20],sigma_D[20][20];
	float M[no_fields-1][no_fields-1];
   FILE* fp=fopen(fname,"r");
   FILE* fp2=fopen(fname1,"w+");
	// printf("Is serial number present in data: \n");
	// int a;
	// scanf("%d",&a);
	int tmp=0;
	float data[MROW][no_fields];
	float classes[no_classes];
	int nos[no_classes];
	memset(nos,0,sizeof(nos));
	int rows=0;
	while(!feof(fp)){
		// if(a){

		// }
		for(int i=0;i<no_fields;i++){
			fscanf(fp,"%f",&data[rows][i]);
			
		}
		rows++;
	}
	if((data[rows-1][0]||data[rows-1][1])==0){ rows-=1;}
	//fclose(fp);
	printf("Data Set reading done\n");
	//printf("Metric No. of data points: %d\n",rows);
	//no_rows[2]=rows;

	for(int i=0;i<rows;i++){
		if(!isPresent(data[i][no_classes-1],classes,tmp)){
			classes[tmp]=data[i][no_classes-1];
			tmp++;
		}
		if(tmp==no_classes) break;
	}
	tmp=0;

	float sum[no_classes][no_fields-1];
	memset(sum,0,sizeof(sum));
	for (int i = 0; i < no_classes; ++i)
	{
		for (int j = 0; j<rows; ++j)
		{
			if(data[j][no_classes-1]==classes[i]){
				nos[i]++;
				for (int k = 0; k < no_fields-1; ++k)
				{
					sum[i][k]+=data[j][k];
				}
			}
		}
	}
   

	// CALCULATING SIGMA S
	float data_subset[MROW][no_fields-1];
	float C[no_fields-1][no_fields-1];
	// float w[no_fields-1];
	// //memset(w,0,sizeof(w));
	// w[0]=-0.5;
	// w[1]=0.25;
	// w[2]=0;
	// w[3]=0;
	int count,full_count=0;
	
	for(int i=0;i<no_classes;i++){
		for(int j=0;j<rows;j++){
			count=0;
			//is class at no_fields-1
			if(classes[i]==data[j][no_fields-1]){
				for(int k=0;k<no_fields-1;k++){
					data_subset[count][k]=data[j][k];
				}
				count++;
			}
		full_count+=count;
		calculate_C(data_subset,C,count);
		//matrix addition
		for (int v=0;v<no_fields-1;v++){
			for(int w=0;w<no_fields-1;w++){
				sum_S[v][w]=sum_S[v][w]+C[v][w];
			}
		}
		//end of matric addition	
		}
		
	}
	//sigma_S=sum_S/full_count;
	//matrix division
		for (int g=0;g<no_fields-1;g++){
			for(int h=0;h<no_fields-1;h++){
				sigma_S[g][h]=sigma_S[g][h]/full_count;
			}
		}
		//end of matric division
	//END OF CALCULATING SIGMA S

	//CALCULATING SIGMA D
	count=0;
	float xi_minus_xj[no_fields-1][no_fields-1],xi_minus_xj_transpose[no_fields-1][no_fields-1];
	float w[no_fields-1];
	//memset(w,0,sizeof(w));
	
	for(int i=0;i<rows;i++){
		//only unique combinations of i and j or not?
		for(int j=0;j<rows;j++){
			
			if(data[i][no_fields-1]!=data[j][no_fields-1]){
			for(int k=0;k<no_fields-1;k++){
				xi_minus_xj[k][k]=data[i][k]-data[j][k];
				transpose(xi_minus_xj,rows,no_fields-1,xi_minus_xj_transpose);
				multiplyMatrices(xi_minus_xj,xi_minus_xj_transpose,C,no_fields-1,no_fields-1,no_fields-1,no_fields-1);
				//matrix addition
				for (int m=0;m<no_fields-1;m++){
					for(int n=0;n<no_fields-1;n++){
						sum_D[m][n]=sum_D[m][n]+C[m][n];
					}
				}
				//end of matric addition
			}
			count++;
		}
		
	}
	}
	//sigma_D=sum_D/count;
	//matrix division
		for (int p=0;p<no_fields-1;p++){
			for(int q=0;q<no_fields-1;q++){
				sigma_S[p][q]=sigma_S[p][q]/full_count;
			}
		}
	//end of matric division
	//END OF CALCULATING SIGMA D
	
	float inverse_S[no_fields-1][no_fields-1],inverse_D[no_fields-1][no_fields-1];
	if(inverse(sigma_S,inverse_S)==true) {
		if(inverse(sigma_D,inverse_D)==true){
			for(int r=0;r<no_fields-1;r++){
				for(int s=0;s<no_fields-1;s++){
					M[r][s]=inverse_S[r][s]-inverse_D[r][s];
				}
			}
		}
	}

	//FILE* fpdat=fopen(fname1,"w+");
	w[0]=-0.5;
	w[1]=0.25;
	w[2]=0;
	w[3]=0;


	float finaldat[rows][2];
	memset(finaldat,0,sizeof(finaldat));
	//printf("%f,%f,%f\n",w[0],w[1],w[2]);
	for(int i=0;i<rows;i++){
		for(int j=0;j<no_fields-1;j++){
			finaldat[i][0]+=(w[j]*data[i][j]);
		}
		finaldat[i][1]=data[i][no_fields-1];
	}

	for(int i=0;i<rows;i++){
		fprintf(fp2,"%f ",finaldat[i][0]);
		fprintf(fp2,"%f\n",finaldat[i][1]);
	}
	fclose(fp);
	fclose(fp2);
	//return;
					
}

void k_nn(char fnames1[3][30]){
	float acc[3];
	float data[3][MROW][2];
	/*float data2[MROW][2];
	float data3[MROW][2];*/

	FILE* fp1;FILE* fp2;FILE* fp3;
	fp1=fopen(fnames1[0],"r");
	fp2=fopen(fnames1[1],"r");
	fp3=fopen(fnames1[2],"r");


	for(int i=0;i<no_rows[0];i++){
		fscanf(fp1,"%f",&data[0][i][0]);
		fscanf(fp1,"%f",&data[0][i][1]);
	}

	for(int i=0;i<no_rows[1];i++){
		fscanf(fp2,"%f",&data[1][i][0]);
		fscanf(fp2,"%f",&data[1][i][1]);
	}

	for(int i=0;i<no_rows[2];i++){
		fscanf(fp3,"%f",&data[2][i][0]);
		fscanf(fp3,"%f",&data[2][i][1]);
	}
	//printf("11111111111111111\n");
	//printf("%f\n",data[1][0][0]);
	

	int k=13;
	float dist[k];
	//memset(dist,-1,sizeof(dist));
	float classes[k];
	//int test_data_indices[160];
	//test_data_indices[0]=4;
	/*for(int i=1;i<160;i++){
		test_data_indices[i]=test_data_indices[i-1]+5;
	}*/
	//int no_test_data;
	float distance;
	float tmp;
	int count[9]={0};
	float correct=0,wrong=0;
	//printf("22222222222222222\n");
	for(int z=0;z<3;z++){
		
		correct=0;
		wrong=0;
		//memset(dist,-1,sizeof(dist));
		//no_test_data=no_rows[z]/5;
		for(int i=4;i<no_rows[z];i=i+5){
			for(int p=0;p<9;p++){
				count[p]=0;
			}
			for(int l=0;l<k;l++){
				dist[l]=-1;
				classes[l]=0;
			}
			for(int j=0;j<no_rows[z];j++){
				if((j+1)%5==0) continue;		//test data point so skip it

				distance=abs(data[z][i][0]-data[z][j][0]);
				//if(z==1&&i==4)printf("%d distance  %f, class %f\n",j,distance,data[z][j][1]);
				if((distance<dist[0])||(dist[0]==-1)){
					//printf("1\n");
					dist[0]=distance;
					classes[0]=data[z][j][1];
					if(z==2) acc[1]<1-acc[1]?acc[1]=1-acc[1]:1;
					for(int l=1;l<k;l++){
						if((dist[l]>dist[l-1])||(dist[l]==-1)){
							tmp=dist[l];
							dist[l]=dist[l-1];
							dist[l-1]=tmp;

							tmp=classes[l];
							classes[l]=classes[l-1];
							classes[l-1]=tmp;
						}
						else break;
					}

				}
				//if(z==1&&i==4)printf("%f %f %f %f %f\n",dist[0],dist[1],dist[2],dist[3],dist[4]);
				//if(z==1&&i==4)printf("%f %f %f %f %f\n",classes[0],classes[1],classes[2],classes[3],classes[4]);	
			}
			//printf("33\n");
			/*for(int l=0;l<k;l++){
				printf("Dist of %d is %f\n",l,dist[l]);
				printf("Class of %d is %f\n",l,classes[l]);
			}*/
			//if(z==1&&i==4)printf("    %d,,,%d\n",count[0],count[1]);
			if(z==1/*||z==0*/)
			for(int l=0;l<k;l++){
				count[l]=0;
			}
			for(int l=0;l<k;l++){
				//if((classes[l]<9)&&(classes[l]>=0))
				(count[(int)(classes[l])])++;
			}
			//printf("44\n");
			//if(z==1&&i==4)printf("%d,,,%d\n",count[0],count[1]);
			int max=0;
			for(int l=1;l<9;l++){
				if(count[(int)max]<=count[l]){
					max=l;
				}
			}
			//if(z==1&&i==4)printf("%d ,, %f\n",max,data[z][i][1]);
			if(max==data[z][i][1]){
				correct++;
			}
			else{
				wrong++;
			}
		}
		acc[z]=correct/(correct+wrong);
		//printf("accuracy for %d is %f\n",z,(correct/(correct+wrong)));
		//printf("%f ,, %f cw\n",correct,wrong);
	}

	
	if(diff_datasets){
		printf("PCA (on Column_2C DataSet)              -- %f\n\nFLD (on ecoli multi-class DataSet)      -- %f\n\nMetric Learning (on transfusion DataSet)-- %f\n\n",acc[0],acc[1],acc[2]);
		printf("Note: Please do not compare Efficiency and Accuracy of three techniques with above accuracies as they are tested against different DataSets\n");
		printf("Please Select Option 2 to compare PCA, FLD and Metric Learning by testing against same DataSet\n\n");
	}
	else
		printf("PCA -- %f\nFLD -- %f\nMetric Learning -- %f\n\n\n",acc[0],acc[1],acc[2]);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
}


int main(int argc,char* argv[]){
	
	no_rows[2]=748;
	//final main code
	
	int no_classes[3]={2,8,2};
	int no_fields[3]={7,8,5};

	char fnames[3][30]={"column_2C.dat","ecoli.data","transfusion.data"};
	char fnames1[3][30]={"pca_res_1","fld_res_1","metric_res_1"};
	char fnames2[3][30]={"pca_res_2","fld_res_2","metric_res_2"};
	int opt;
	while(1){
	printf("If you want PCA FLD and Metric Learning to run on 3 different DataSets, Please Select Option 1\n");
	printf("If you want all three to run on same DataSet (to compare their accuracies) Please Select Option 2\n");
	printf("If you want to quit program Select Option 0\n");
	printf("Please Select Option 0 or 1 or 2: ");
	scanf("%d",&opt);
	if(opt==1){

		diff_datasets=true;
		printf("\n\nStarting PCA\n\n");
		pca(no_classes[0],no_fields[0],fnames[0],fnames1[0]);
		printf("\nPCA Done!!\n\n");
		printf("----------------------------------------------------");
		printf("\n\nStarting FLD\n\n");
		fld(no_classes[1],no_fields[1],fnames[1],fnames1[1]);
		printf("\nFLD Done!!\n\n");
		printf("----------------------------------------------------");
		printf("\n\nStarting Metric Learning\n\n");
		metric_learning(no_classes[2],no_fields[2],fnames[2],fnames1[2]);
		printf("\nMetric Learning Done!!\n\n");
		printf("----------------------------------------------------\n");

		printf("Accuracies of each models on different DataSets:\n");
		k_nn(fnames1);
	}
	if(opt==2){

		diff_datasets=false;
		printf("\n\nStarting PCA\n\n");
		pca(no_classes[2],no_fields[2],fnames[2],fnames2[0]);
		printf("\nPCA Done!!\n\n");
		printf("----------------------------------------------------");
		printf("\n\nStarting FLD\n\n");
		fld(no_classes[2],no_fields[2],fnames[2],fnames2[1]);
		printf("\nFLD Done!!\n\n");
		printf("----------------------------------------------------");
		printf("\n\nStarting Metric Learning\n\n");
		metric_learning(no_classes[2],no_fields[2],fnames[2],fnames2[2]);
		printf("\nMetric Learning Done!!\n\n");
		printf("----------------------------------------------------\n");

		printf("Accuracies on transfusion DataSet:\n");
		k_nn(fnames2);
	}
	if(opt==0){break;}
	else{
		continue;
	}
	}
	

	
	
