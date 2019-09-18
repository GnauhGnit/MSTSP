# Multi-solution Traveling Salesman Problem (MSTSP) 


## 1. Introduction

Multi-solution Traveling Salesman Problem (MSTSP) is essentially a TSP, but the one with multiple optimal solutions. This benchmark included 25 MSTSPs. The number of cities ranges from 9 to 66, and the number of optimal solutions ranges from 4 to 196. The algorithms used to solve MSTSPs are required to provide a solution set. These candidate solutions are further evaluated with two indicators, i.e., Fbeta and Diversity Indicator (DI). Fbeta measures the solution quality, and DI measures the solution diversity. 

## 2. Documents

In this directory, there are two folders. The benchmark_MSTSP folder contains MSTSP instances and their optimal solutions. The demo folder provides a demo implemented in MATLAB to show how to evaluate the obtained solutions. The details are listed in the following. 

The benchmark_MSTSP folder includes the following files: <br>
**<instance_name>.tsp** <br>
**<instance_name>.solution** <br>
in which: <br>
	**<instance_name>:** <br>
is the name of MSTSP instance corresponding to the class, the index and the number of city of the instance. For example, the name of the first instance, which is a simple MSTSP with 9 cities, names **simple1_9.tsp** <br>
**<instance_name>.tsp** <br>
		The file contains the two-dimensional coordinates of the cities in **simple1_9**. Every row has two numbers, which are the abscissa and the ordinate. The number of rows is exactly the number of cities. <br>
**<instance_name>.solution** <br>
		The file contains all the optimal solutions of <instance_name>. Every row lists one optimal solution. The first number is the tour length of this tour, and the remaining numbers construct the permutation of this optimal solution. The number of rows represents the number of optimal solutions of <instance_name>. <br>

The demo folder includes the following MATLAB files:<br>
**measure_share_dist.m** <br>
	The file is a MATLAB code to compute the share distance between two solutions. <br>
**MSTSP_measure.m** <br>
		The file is a MALAB code to input obtained solutions and <instance_name> index, and then return their Fbeta and DI. <br>
**MSTSP_demo.m** <br>
		The file is a MATLAB code to demonstrate how to evaluate benchmark suite of MSTSP. 

## 3. How to use
	
Run the **MSTSP_demo.m**, and you will see the evaluation results for randomly generated solutions. 

------

\* For more information please refer to the Technical Report of Multi-solution Traveling Salesman Problem [1]. Besides, a niching memetic algorithm is proposed to solve the MSTSPs [2].  Any comments and suggestions are welcome.

[1] T. Huang, Y.-J. Gong and J. Zhang, “Seeking Multiple Solutions of Combinatorial Optimization Problems: A Proof of Principle Study,” in 2018 IEEE Symposium Series on Computational Intelligence (SSCI), 2018.
[2] T. Huang, Y. Gong, S. Kwong, H. Wang and J. Zhang, “A Niching Memetic Algorithm for Multi-Solution Traveling Salesman Problem,” IEEE Transactions on Evolutionary Computation. DOI: 10.1109/TEVC.2019.2936440.