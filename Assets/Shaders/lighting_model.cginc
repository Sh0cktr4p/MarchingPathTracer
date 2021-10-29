
interface LightingModel {
    fixed3 eval(float3 normal, float3 in_dir, float3 out_dir);
};


class TorranceSparrow : LightingModel {
    float _roughness;
    float3 _material;

    fixed3 eval(float3 normal, float3 in_dir, float3 out_dir) {
        float3 l = in_dir;
        float3 v = -out_dir;
        float3 h = (l + v) / 2;
        h = normalize(h);
        float3 n = normal;

        _material = float3(1.00, 0.71, 0.29); // Gold
        //float3 f0 = float3(0.91, 0.92, 0.92); // Aluminium
        //float3 f0 = float3(0.95, 0.64, 0.54); // Copper
        float3 f = _material + (1 - _material) * pow(1 - dot(n, l), 5);

        float betha = acos(dot(h, n));
        float d = exp(-pow(tan(betha), 2) / pow(roughness, 2)) / (3.1416926 * pow(roughness, 2) * pow(max(dot(h, n), 0.0001), 4));
        float g = min(1, min(2 * dot(n, h) * dot(n, v) / dot(v, h), 2 * dot(n, h) * dot(n, l) / dot(v, h)));
        if (isnan(d)) {
            return float4(1, 0, 0, 1);
        }
        float3 result = saturate(f * d * g / (3.1415926 * dot(n, l) * dot(n, v))) * dot(n, l);
        return result;
    }
};


fixed4 phong(float3 normal, float3 in_dir, float3 out_dir, float roughness) {
    float3 n = normal;
    float3 l = in_dir;
    float3 v = out_dir;

    float diff = saturate(dot(n, l));
    float spec = saturate(pow(dot(v, -reflect(v, n)), 5));

    float3 result = saturate(0.1 + diff * roughness + spec * (1 - roughness));
    return fixed4(result, 1);
}

fixed4 blinn_phong(float3 normal, float3 in_dir, float3 out_dir, float roughness) {
    float3 n = normal;
    float3 l = in_dir;
    float3 v = -out_dir;
    float3 h = (v + l) / 2;
    h = normalize(h);
    float hn = dot(h, n);

    float ref = hn <= 0 ? 0 : saturate(pow(hn, roughness));
    float3 result = saturate(0.01 + ref);
    return fixed4(result, 1);
}

fixed4 torrance_sparrow(float3 normal, float3 in_dir, float3 out_dir, float roughness) {
    float3 l = in_dir;
    float3 v = -out_dir;
    float3 h = (l + v) / 2;
    h = normalize(h);
    float3 n = normal;

    float3 f0 = float3(1.00, 0.71, 0.29); // Gold
    //float3 f0 = float3(0.91, 0.92, 0.92); // Aluminium
    //float3 f0 = float3(0.95, 0.64, 0.54); // Copper
    float3 f = f0 + (1 - f0) * pow(1 - dot(n, l), 5);

    float betha = acos(dot(h, n));
    float d = exp(-pow(tan(betha), 2) / pow(roughness, 2)) / (3.1416926 * pow(roughness, 2) * pow(max(dot(h, n), 0.0001), 4));
    float g = min(1, min(2 * dot(n, h) * dot(n, v) / dot(v, h), 2 * dot(n, h) * dot(n, l) / dot(v, h)));
    if (isnan(d)) {
        return float4(1, 0, 0, 1);
    }
    float3 result = saturate(f * d * g / (3.1415926 * dot(n, l) * dot(n, v))) * dot(n, l);
    return float4(result, 1);
}