
#define GOLD float3(1.00, 0.71, 0.29)
#define ALUMINIUM float3(0.91, 0.92, 0.92)
#define COPPER float3(0.95, 0.64, 0.54)


interface LightingModel {
    fixed3 eval(float3 normal, float3 in_dir, float3 out_dir);
};


class TorranceSparrow : LightingModel {
    float _roughness;
    float3 _material;

    fixed3 eval(float3 normal, float3 in_dir, float3 out_dir) {
        float3 l = in_dir;
        float3 v = out_dir;
        float3 h = (l + v) / 2;
        h = normalize(h);
        float3 n = normal;

        float3 f = _material + (1 - _material) * pow(1 + dot(n, l), 5);

        float betha = acos(dot(h, n));
        float d = exp(-pow(tan(betha), 2) / pow(_roughness, 2)) / (3.1416926 * pow(_roughness, 2) * pow(dot(h, n)+ 0.0001, 4));
        float g = min(1, min(2 * dot(n, h) * dot(n, v) / dot(v, h), 2 * dot(n, h) * dot(n, l) / dot(v, h)));
        if (isinf(d)) {
            return float4(1, 0, 0, 1);
        }
        float3 result = saturate(f * g * d/ (3.1415926 * dot(n, l) * dot(n, v)));
        return result * dot(n, l);
    }
};

TorranceSparrow torrance_sparrow(float roughness, float3 material) {
    TorranceSparrow ts = {
        roughness,
        material
    };

    return ts;
}

class Lambert : LightingModel {
    float3 _material;

    fixed3 eval(float3 normal, float3 in_dir, float3 out_dir) {
        return _material * dot(-out_dir, normal);
    }
};

Lambert lambert(float3 material) {
    Lambert l = {
        material
    };

    return l;
}


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
