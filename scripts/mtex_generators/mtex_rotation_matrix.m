function matrixValue = mtex_rotation_matrix(rotationLike)
%MTEX_ROTATION_MATRIX Convert an MTEX orientation/rotation to a 3x3 matrix.
try
    matrixValue = matrix(rotationLike);
catch
    matrixValue = double(rotationLike);
end
end

