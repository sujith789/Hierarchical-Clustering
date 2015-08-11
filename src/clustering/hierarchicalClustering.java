package clustering;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;


public class hierarchicalClustering {
	
	public static void main(String args[]) throws Exception
	{	
		String filename=null;
		int mGenes=0;
		System.out.println("Enter the file");
		BufferedReader reader=new BufferedReader(new InputStreamReader(System.in));
         filename = reader.readLine();
		
		BufferedReader readfile = new BufferedReader(new FileReader(filename));
		while(readfile.readLine()!=null)	
		{
			mGenes++;
		}
		//System.out.println("Total no of transactions ="+mGenes);
		int groundtruth[] = new int[mGenes];
		
		BufferedReader readfile1 = new BufferedReader(new FileReader(filename));
		String read;
		read=readfile1.readLine();
		String tokens[]=read.split("\\t");
		int columns =tokens.length-2;
		
		//System.out.println(" Number of columns = "+columns);
		
		double data_matrix[][] = new double[mGenes][columns];
		
		
		BufferedReader readfile2 = new BufferedReader(new FileReader(filename));
		String read1;
		int row1=0;
		int column1=0;
		while ((read1 = readfile2.readLine()) != null) 
		{
			//read1=readfile2.readLine();
			String tokens1[]=read1.split("\\t");
			column1=0;
			for(int i=0;i<=1;i++)
			{
				String geneid=tokens1[0];
				String truth=tokens1[1];
				int gid=Integer.parseInt(geneid);
				int gtruth=Integer.parseInt(truth);
				groundtruth[row1]= gtruth;
			}
			for(int i=2;i<tokens1.length;i++)
			{
				double value=Double.parseDouble(tokens1[i]);
				//System.out.println(value);
				data_matrix[row1][column1]=value;
				column1++;
				
			}
			row1++;
		}
		
		double distance_matrix[][] = new double[mGenes][mGenes];
		double new_distance_matrix[][]= new double[mGenes][mGenes];
		
		System.out.println();
		for(int i=0;i<row1;i++)
		{
			for(int j=0;j<row1;j++)
			{
				if(i==j)
				{
					distance_matrix[i][j]=0.0;
					new_distance_matrix[i][j]=0.0;
				}
				else
				{	double sum=0;
				    double temp=0;
					for(int k=0;k<column1;k++)
					{  
						temp=data_matrix[i][k]-data_matrix[j][k];
						double sq_temp=temp*temp;
						sum=sum+sq_temp;
					}
					double root_temp=Math.sqrt(sum);
					distance_matrix[i][j]=root_temp;
					new_distance_matrix[i][j] = root_temp;
					
				}
			}
			//System.out.println();
		}//end of for
		
		//printing the distance matrix;
		
		/*System.out.println(" The distance matrix is as follows..");
		System.out.println();
		for(int i=0;i<row1;i++)
		{
			for(int j=0;j<row1;j++)
			{
				System.out.print(distance_matrix[i][j] + "  ");
			}
			
		}*/
		
		//new_distance_matrix = distance_matrix;
		
		int count_values=row1;
		int clusternum=1;
		
		Map<Integer,Integer> clustermap = new TreeMap<Integer, Integer>();
		for(int m=0;m<mGenes;m++)
		{
			clustermap.put(m, m);
		}
		int ground_value = num_elements(groundtruth);
		System.out.println(" Ground value = "+ground_value);
		
		while(count_values>ground_value)
		{
		double smallest = 100;
		int srow=-1,scol=-1;
		
		
		for(int i=0;i<row1;i++)
		{
			for(int j=0;j<row1;j++)
			{
				if((smallest>distance_matrix[i][j]) && (distance_matrix[i][j]!=0))
						{
					     smallest = distance_matrix[i][j];
					     srow=i;
					     scol=j;
						}
			}
			
		}
		
		
		//for finding min distance
		for(int k=0;k<row1;k++)
		{
			if(distance_matrix[srow][k]>distance_matrix[scol][k])
			{
				distance_matrix[srow][k]=distance_matrix[scol][k];
				distance_matrix[k][srow] = distance_matrix[scol][k];
			}
		}
		
		//for finding the group average
		/*for(int k=0;k<row1;k++)
		{   if((k!=srow) ||(k!=scol))
		{
			distance_matrix[srow][k]= (distance_matrix[srow][k]+distance_matrix[scol][k])/2;
			distance_matrix[k][srow] =distance_matrix[srow][k];
		}
		}*/
		
		
		
		//updating all the values in the row of the column which corresponds to the smallest value
		for(int k=0;k<row1;k++)
		{
			distance_matrix[scol][k]=0;
			distance_matrix[k][scol]=0;
		}
		
		if(clustermap.containsKey(srow)||clustermap.containsKey(scol))
		{
			
			if((clustermap.containsKey(srow))&&(clustermap.containsKey(scol)))
			{
				int getnum = clustermap.get(srow);
				int getcol = clustermap.get(scol);
				clustermap.put(scol,getnum);
				for(Map.Entry<Integer, Integer> entry: clustermap.entrySet())
				{
					if(entry.getValue()==getcol)
					{
						clustermap.put(entry.getKey(), getnum);
					}
				}
			}
			else if(clustermap.containsKey(srow))
			{
				int getnum = clustermap.get(srow);
				clustermap.put(scol, getnum);
			}
			else if(clustermap.containsKey(scol))
			{  
				int getnum = clustermap.get(scol);
				clustermap.put(srow, getnum);
			}
		}
		else
		{
			clustermap.put(srow, clusternum);
			clustermap.put(scol, clusternum);
			clusternum++;
		}
		
		count_values--;
		}//end of while
		
		int foundtruth[] = new int[mGenes];
		System.out.println(" Gene Id and corresponding cluster number is as follows -");
		
		for(Map.Entry<Integer, Integer> entry: clustermap.entrySet())
		{
			//if(entry.getValue()!=71)
			foundtruth[entry.getKey()]=entry.getValue();
			//System.out.println((entry.getKey()+1) +"\t" + (entry.getValue()+1));
			System.out.println(" Gene Id = "+(entry.getKey()+1) + " Cluster number = "+(entry.getValue()+1));
		}
		
		
		//Incident_matrix calculation
		int incident_matrixP[][] = new int[mGenes][mGenes];
		int incident_matrixC[][] = new int[mGenes][mGenes];
		int ss=0,dd=0,sd=0,ds=0;
		for(int i=0;i<mGenes;i++)
		{
			for(int j=0; j<mGenes; j++)
			{
				if(groundtruth[i]==groundtruth[j])
				{
					incident_matrixP[i][j] =1;
				}
				
				else
				{
					incident_matrixP[i][j]=0;
				}
				
				if(foundtruth[i]==foundtruth[j])
				{
					incident_matrixC[i][j]=1;
				}
				else
				{
					incident_matrixC[i][j]=0;
				}
				
			}
		}//end of for
		
	for(int i=0;i<mGenes;i++)	
	{
		for(int j=0;j<mGenes;j++)
		{
			if(incident_matrixP[i][j]==incident_matrixC[i][j])
			{
				if(incident_matrixP[i][j]==1)
				{
					ss++;
				}
				else
				{
					dd++;
				}
			}
			else
			{
				if(incident_matrixP[i][j]==1)
				{
					ds++;
				}
				else
				{
					sd++;
				}
			}
		}
	}
		
	System.out.println("SS = "+ss);
	System.out.println("DD = "+dd);
	System.out.println("SD = "+sd);
	System.out.println("DS = "+ds);
	
	double rand_index;
	double jaccard_coef;
	double ss1=(double)ss;
	double dd1 = (double)dd;
	double sd1= (double)sd;
	double ds1 = (double)ds;
	
	rand_index =(ss1+dd1)/(ss1+dd1+sd1+ds1);
	jaccard_coef=(ss1)/(ss1+sd1+ds1);
	
	System.out.println(" Rand index = "+rand_index);
	System.out.println(" Jacard Coefficient = "+jaccard_coef);
	
	//internal_index calculation
	double cincident_matrix[][] = new double[mGenes][mGenes];
	for(int i=0;i<mGenes;i++)
	{
		for(int j=0; j<mGenes;j++)
		{
			cincident_matrix[i][j] =(double)incident_matrixC[i][j];
		}
	}
	
	double dmean = mean(new_distance_matrix,mGenes);
	double cmean = mean(cincident_matrix,mGenes);
	double nvar=numerical_variance(new_distance_matrix, cincident_matrix, mGenes);
	//System.out.println(" Nvar = "+nvar);
	double dvar=Math.sqrt(variance(new_distance_matrix,mGenes));
	//System.out.println(" Dvar = "+dvar);
	double cvar = Math.sqrt(variance(cincident_matrix, mGenes));
	//System.out.println(" Cvar = "+cvar);
	double correlation = Math.abs(nvar/(dvar*cvar));
	System.out.println(" Correlation of incident matrix and distance matrix = "+correlation);
	}
	
	
	public static double mean(double arr[][], int num)
	{
		double sum=0,mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
			}
		}
		mean = sum/(num*num);
			
		return mean;
	}
	
	public static double variance(double arr[][], int num)
	{
		double sum=0,mean=0, var=0,sum_mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
			}
		}
		mean = sum/(num*num);
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum_mean= sum_mean + ((arr[i][j]-mean)*(arr[i][j]-mean));
			}
		}
		return sum_mean;
	}
	
	public static double numerical_variance(double arr[][],double arr1[][],int num)
	{
		double sum=0,sum1=0,mean=0,mean1=0,sum_mean=0;
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum=sum+arr[i][j];
				sum1=sum1+arr1[i][j];
			}
		}
		mean = sum/(num*num);
		mean1=sum1/(num*num);
		for(int i=0;i<num;i++)
		{
			for( int j=0;j<num;j++)
			{
				sum_mean= sum_mean + ((arr[i][j]-mean)*(arr1[i][j]-mean1));
			}
		}
		return sum_mean;
	}
	
	public static int num_elements(int[] arr)
	{
		Set<Integer> newset = new HashSet<Integer>();
		for(int element: arr)
		{
			newset.add(element);
		}
		if(newset.contains(-1))
		{
			return newset.size()-1;
		}
		else
		{
			return newset.size();
		}
		
	}
	
}
