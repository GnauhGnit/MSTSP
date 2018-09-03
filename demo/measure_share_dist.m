function [ share_dist ] = measure_share_dist( vec1, vec2 )
%   return the share dist
% vec1 = [1 2 3 4 5]
% vec2 = [5 4 3 2 1]
len = size(vec1, 2);
share_dist = 0;
new_vec1 = vec1'* ones(1, 2);
new_vec2 = vec2'* ones(1, 2);

similar_matrix_1 = zeros(len, len);
similar_matrix_2 = zeros(len, len);

new_vec1 = reshape(new_vec1', length(new_vec1(:)), 1)';
new_vec1 = [new_vec1(2:end) new_vec1(1)];
new_vec1 = reshape(new_vec1', 2, length(new_vec1(:))/2 )';
similar_matrix_1(sub2ind(size(similar_matrix_1), new_vec1(:, 1),  new_vec1(:, 2)))=1;
similar_matrix_1(sub2ind(size(similar_matrix_1), new_vec1(:, 2),  new_vec1(:, 1)))=1;

new_vec2 = reshape(new_vec2', length(new_vec2(:)), 1)';
new_vec2 = [new_vec2(2:end) new_vec2(1)];
new_vec2 = reshape(new_vec2', 2, length(new_vec2(:))/2 )';
similar_matrix_2(sub2ind(size(similar_matrix_2), new_vec2(:, 1),  new_vec2(:, 2)))=1;
similar_matrix_2(sub2ind(size(similar_matrix_2), new_vec2(:, 2),  new_vec2(:, 1)))=1;

share_dist = sum(sum(similar_matrix_1&similar_matrix_2))/2;
end

