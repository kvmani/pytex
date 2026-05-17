function q = mtex_quaternion_wxyz(rotationLike)
%MTEX_QUATERNION_WXYZ Convert an MTEX orientation/rotation to [w x y z].
quat = quaternion(rotationLike);
q = double([quat.a, quat.b, quat.c, quat.d]);
if q(1) < 0
    q = -q;
end
q = q ./ norm(q);
end

