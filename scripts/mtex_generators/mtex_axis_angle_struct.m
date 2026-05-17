function payload = mtex_axis_angle_struct(rotationLike)
%MTEX_AXIS_ANGLE_STRUCT Return axis-angle representation.
try
    axisValue = axis(rotationLike);
    angleValue = angle(rotationLike);
catch
    rot = rotation(rotationLike);
    axisValue = axis(rot);
    angleValue = angle(rot);
end
payload = struct( ...
    "axis", mtex_vector_xyz(axisValue), ...
    "angle_deg", double(angleValue / degree));
end

