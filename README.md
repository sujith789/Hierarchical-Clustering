# Hierarchical-Clustering

Problem Description
--------------------------
To implement hierarchical clustering to find clusters of genes, which exhibit similar expression profiles.

How to Run
------------------------
1) Run HierarchicalClustering.java
2) Enter the desired file name - "cho.txt" or "iyer.txt". 
3) Output is displayed.

Flow of algorithm
----------------------------
•	Initially, all the points are in one cluster.
•	Calculate the number of distinct clusters from ground truth value, and call it say ‘V’
•	Let the total number of items be ‘N’.
•	Create a new cluster map.
•	Repeat the following steps, N- V times
•	Choose the minimum value from the distance matrix. Let the minimum value be in row ‘i’ and column’j’.
•	Check any other values in the cluster map have same value as ‘i’. If so update all those values including ‘i’ and put in the new cluster.
•	Do the same thing for column ‘j’ as well. 
•	This means that all the values in new cluster which have same values as ‘i’ or ‘j’, together with ‘i’ as well as ‘j’ goes to the same cluster.
•	After we have done the above steps for N-V times, the resultant new cluster will contains only that many clusters.
•	The values assigned to new cluster is the final values

