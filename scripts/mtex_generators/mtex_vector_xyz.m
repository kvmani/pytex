function xyz = mtex_vector_xyz(vectorValue)
%MTEX_VECTOR_XYZ Convert an MTEX vector3d to a numeric row.
xyz = double([vectorValue.x, vectorValue.y, vectorValue.z]);
end

