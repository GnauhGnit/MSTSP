% ******************************************************************************
% * Version: 1.0
% * Last modified on: 28 May, 2018
% * Developers: Ting Huang
% *      email: cshting@mail.scut.edu.cn
% * ****************************************************************************

% Demonstration file on how to evaluate the benchmark suite of MSTSP

clear;
clc;

CITY_SIZE = [9	10	10	11	12	12	10	12	10	10	10	15	28	34	22	33	35	39	42	45	48	55	59	60	66];
Max_FEs = [60000*ones(1,12) 1200000*ones(1,13)]; 
NP = 100;

for index = 1:25
    alg_solution = zeros(NP, CITY_SIZE(index));
    for j = 1 :NP
        alg_solution(j, :) = randperm(CITY_SIZE(index));
    end
    
    [Fbeta, DI] = MSTSP_measure(alg_solution, index);
    fprintf('%dth MSTSP: F_beta = %f, DI = %f\n', index,Fbeta, DI);
end


