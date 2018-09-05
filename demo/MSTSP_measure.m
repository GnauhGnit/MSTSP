function [Fbeta, DI] = MSTSP_measure(alg_solution, index)
MSTSP_NAME = { 'simple1_9', 'simple2_10', 'simple3_10', 'simple4_11', 'simple5_12', 'simple6_12', ...
    'geometry1_10', 'geometry2_12', 'geometry3_10', 'geometry4_10', 'geometry5_10', 'geometry6_15', ...
    'composite1_28','composite2_34','composite3_22','composite4_33','composite5_35','composite6_39','composite7_42','composite8_45', ...
    'composite9_48','composite10_55','composite11_59','composite12_60','composite13_66'};
% MSTSP_OPT_LEN = [680	1265	832	803	754	845	130	1344	72	72	78	130	3055	3575	9455	8761	9061	23763	14408	10973	6767	10442	24451	9614	9521];
% MSTSP_OPT_SIZE = [3	4	13	4	2	4	56	110	4	4	14	196	70	16	72	64	10	20	20	20	4	9	10	36	26];
MSTSP_BASEPATH = '../benchmark_MSTSP/';
beta2 = 0.3;

if isempty(alg_solution)
    fprintf('The number of the solutions is empty.\n');
    return
elseif index < 1 || index > 25
    fprintf('The index is out of range(1-25).\n');
    return
end
mstsp_cities_cardinate = load(strcat(MSTSP_BASEPATH, char(MSTSP_NAME(index)), '.tsp')); 
mstsp_solution = load(strcat(MSTSP_BASEPATH,  char(MSTSP_NAME(index)), '.solution')) + 1;
mstsp_solution = mstsp_solution(:, 1:end-1);

% calculate the length of tours
alg_solution = [zeros(size(alg_solution, 1), 1)  alg_solution];
for i = 1: size(alg_solution, 1)
    pathtour = mstsp_cities_cardinate([alg_solution(i, 2:end) alg_solution(i, 2)], :);
    alg_solution(i, 1) = sum(round(sqrt(sum((pathtour(1:end - 1, :) - pathtour(2:end, :)).^2, 2))));
end

city_num = size(mstsp_solution, 2) - 1;
share_dist = zeros(size(mstsp_solution, 1), size(alg_solution, 1));
for i = 1:size(mstsp_solution, 1)
    for j = 1: size(alg_solution, 1)
        share_dist(i, j) = measure_share_dist(alg_solution(j, 2:end), mstsp_solution(i, 2:end));
    end
end

flag_mstsp = max(share_dist, [], 2) == repmat(city_num,size(mstsp_solution, 1), 1 );
flag_alg = max(share_dist, [], 1)' == repmat(city_num,size(alg_solution, 1), 1 );
[nearest_num]  = max(share_dist, [], 2);


DI = mean(nearest_num) / city_num;
TP =  sum(flag_alg);
FP =  size(flag_alg, 1) - TP;
FN = size(flag_mstsp, 1) - sum(flag_mstsp);
P = TP / (TP + FP);
R = TP / (TP + FN);
Fbeta = (1+beta2)*P*R /((beta2)*P + R);
if isnan(Fbeta)
    Fbeta = 0;
end

end