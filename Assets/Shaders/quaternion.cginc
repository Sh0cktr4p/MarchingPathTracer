typedef float4 quaternion;

quaternion quat(float3 axis, float angle){
    return quaternion(axis * sin(angle / 2), cos(angle / 2));
}

quaternion conj(quaternion q){
    return quaternion(-q.xyz, q.w);
}

quaternion qmul(quaternion a, quaternion b){
    return quaternion(a.xyz * b.w + b.xyz * a.w + cross(a.xyz, b.xyz), a.w * b.w - dot(a.xyz, b.xyz));
}

float3 qrot(quaternion q, float3 v){
    return qmul(qmul(q, quaternion(v, 0)), conj(q)).xyz;
}