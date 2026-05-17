function array = mtex_vector_array_xyz(vectorValues)
%MTEX_VECTOR_ARRAY_XYZ Convert an MTEX vector3d array to n-by-3 numeric rows.
count = numel(vectorValues);
array = zeros(count, 3);
for idx = 1:count
    array(idx, :) = mtex_vector_xyz(vectorValues(idx));
end
end

