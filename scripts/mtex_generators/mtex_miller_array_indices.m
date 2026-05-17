function array = mtex_miller_array_indices(millerValues)
%MTEX_MILLER_ARRAY_INDICES Convert an MTEX Miller array to numeric rows.
count = numel(millerValues);
array = zeros(count, 3);
for idx = 1:count
    array(idx, :) = mtex_miller_indices(millerValues(idx));
end
end

