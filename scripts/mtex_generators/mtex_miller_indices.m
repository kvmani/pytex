function indices = mtex_miller_indices(millerValue)
%MTEX_MILLER_INDICES Best-effort conversion of an MTEX Miller value to hkl/uvw triplet.
try
    indices = double([millerValue.h, millerValue.k, millerValue.l]);
catch
    indices = double(millerValue);
end
end

