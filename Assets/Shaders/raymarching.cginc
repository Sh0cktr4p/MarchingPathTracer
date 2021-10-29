#include "sdf.cginc"
#include "lighting_model.cginc"

struct Ray {
    float3 _pos;
    float3 _dir;
};

Ray to_ray(float3 pos, float3 dir) {
    Ray ray = {
        pos,
        dir
    };

    return ray;
}


struct MarchResult{
    float3 _pos;
    float _min_dist;
    int _steps;
};

MarchResult to_result(float3 pos, float min_dist, int steps){
    MarchResult result = {
        pos, 
        min_dist, 
        steps
    };

    return result;
}


class RayMarcher {
    float _eps;
    float _max_dist;
    int _max_steps;

    bool hit_something(MarchResult result) {
        return result._min_dist < _eps;
    }

    bool hit_sdf(SDF sdf, MarchResult result) {
        return sdf.eval(result._pos) < _eps;
    }

    MarchResult march(SDF sdf, Ray ray) {
        float dist = _eps;
        float min_dist = 1.0; // TODO check if higher values are better

        float3 curr_pos = ray._pos + ray._dir * dist;

        for (int i = 0; i < _max_steps; i++) {
            float new_dist = sdf.eval(curr_pos); // TODO here stood "* 0.4"
            dist += new_dist;
            //curr_pos = ray._pos + ray._dir * dist;
            curr_pos += ray._dir * new_dist;

            min_dist = min(min_dist, 3 * new_dist / dist);

            if (new_dist < _eps) {
                // Hit something
                return to_result(curr_pos, 0, i);
            }

            if (dist > _max_dist) {
                // Left the scene
                return to_result(ray._pos, min_dist, i);
            }
        }

        // Out of steps, should ideally be the exception
        return to_result(ray._pos, min_dist, max_steps);
    }

    float get_dist(SDF sdf, Ray ray) {
        float dist = _eps;
        float min_dist = 1.0;

        float3 curr_pos = ray._pos + ray._dir * dist;

        for (int i = 0; i < _max_steps; i++) {
            float new_dist = sdf.eval(curr_pos);
            dist += new_dist;
            //curr_pos = pos + dir * dist;
            curr_pos += ray._dir * new_dist;

            if (new_dist < _eps) {
                return dist;
            }

            if (dist > _max_dist) {
                return dist;
            }
        }

        return dist;
    }

    float3 calc_normal(SDF sdf, float3 pos) {
        float2 e2 = float2(_eps, 0.0);

        float3 nor = float3(
            sdf.eval(pos + e2.xyy) - sdf.eval(pos - e2.xyy),
            sdf.eval(pos + e2.yxy) - sdf.eval(pos - e2.yxy),
            sdf.eval(pos + e2.yyx) - sdf.eval(pos - e2.yyx)
        );

        return normalize(nor);
    }
};

RayMarcher new_ray_marcher(float eps, float max_dist, float max_steps) {
    RayMarcher ray_marcher = {
        eps,
        max_dist,
        max_steps
    };

    return ray_marcher;
}



// Deprecated stuff //


float softshadow(SDF sdf, float3 ro, float3 rd, float mint, float maxt, float k){
    float res = 1.0;
    float ph = 1e20;
    for(float t = mint; t < maxt;){
        float h = sdf.eval(ro + rd * t);
        if(h < 0.001){
            return 0.0;
        }
        float y = h * h / (2.0 * ph);
        float d = sqrt(h * h - y * y);
        res = min(res, k * d / max(0.0, t - y));
        ph = h;
        t += h;
    }
    return res;
}

float ambient_occlusion(SDF sdf, float3 pos, float3 nor){
	float occ = 0.0;
    float sca = 1.0;
    for(int i=0; i<5; i++){
        float h = 0.001 + 0.15*float(i)/4.0;
        float d = sdf.eval( pos + h*nor );
        occ += (h-d)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 1.5*occ, 0.0, 1.0 );    
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


fixed4 local_illumination(SDF sdf, float3 pos, float3 dir, float3 col, float3 ild, float3 lco, float ambient, float intensity, 
        float spec_pow, int shadow_pow, float a_occ, float step_count, float epsilon, float max_dist, int max_steps, float dist){
    float3 n = calcNormal(sdf, pos, epsilon);
    float3 l = normalize(ild);
    float3 v = normalize(dir);
    l = -v;
    float3 h = normalize(l - v);

    float nl = saturate(dot(n, l));

    float3 diff = saturate(max(ambient, nl) * col * lco);
    float3 blin = saturate(pow(saturate(dot(h, n)), spec_pow)) * lco;

    float3 result = diff + blin * 0.1;

    //float shadow_value = softshadow(sdf, pos, l, 2 * epsilon, max_dist, shadow_pow);
    //float shadow_value = march(sdf, pos + n * 2 * epsilon, l, epsilon, max_dist, max_steps).min_dist;
    //shadow_value = pow(saturate(shadow_value), shadow_pow);

    //result = max(result * softshadow(sdf, pos, l, 0.01, dist, shadow_pow), ambient);

    // Ambient occlusion
    //result *= (1 - saturate(step_count / a_occ));
    //result *= ambient_occlusion(sdf, pos, n);

    float glow_pulse = 0*max(0, sin(_Time.z + 3.1415926*pow(1 - dist / max_dist, 2)) - 0.9) * 10;

    //return phong(n, l, v, 0.);
    //return blinn_phong(n, l, v, spec_pow);
    return torrance_sparrow(n, l, v, spec_pow);
    return fixed4(result * (pow(min(intensity / dist, 1.4), 2) + glow_pulse), 1);
}



fixed4 blinn_phong_raymarch(SDF sdf, float3 pos, float3 dir, float3 inv_light_dir, float4 back_col, float3 obj_col, float3 light_col, 
        float ambient, float intensity, float spec_pow, int shadow_pow, float a_occ, float epsilon, float max_dist, int max_steps){
    MarchResult march_result = march(sdf, pos, dir, epsilon, max_dist, max_steps);

    fixed4 final_col = back_col;
    final_col = 0;

    float dist = length(march_result.pos - pos);

    if(dist > epsilon){
        final_col = local_illumination(sdf, march_result.pos, dir, obj_col, inv_light_dir, light_col, ambient, intensity, spec_pow, shadow_pow, a_occ, march_result.steps, epsilon, max_dist, max_steps, dist);
    }

    return final_col;
}

fixed4 menger_lod_march(float3 pos, float3 dir, float3 inv_light_dir, float4 back_col, float3 obj_col, float3 light_col, 
        float ambient, float intensity, int spec_pow, int shadow_pow, float a_occ, float epsilon, float max_dist, int max_steps,
        int menger_lvl, float menger_scl, float fold_mnr, float fold_mxr, float rep_scl, float lim){
    
    SDF basic_sdf = repeat(sphere_fold(menger_sponge(menger_lvl, menger_scl), fold_mnr, fold_mxr), rep_scl);
    float dist = april(basic_sdf, pos, dir, epsilon, max_dist, max_steps);
        
    fixed4 final_col = back_col;

    if(dist > epsilon && dist < max_dist){
        int final_lvl = min(menger_lvl, log10(menger_scl / (lim * dist)) / log10(3));
        SDF final_sdf = repeat(sphere_fold(menger_sponge(final_lvl, menger_scl), fold_mnr, fold_mxr), rep_scl);

        MarchResult final_result = march(final_sdf, pos, dir, epsilon, max_dist, max_steps);
    
        dist = length(final_result.pos - pos);

        //obj_col = lerp(float3(1, 0, 0), float3(0, 0, 1), exp(-2.71*dist / max_dist));

        if(dist > epsilon){
            final_col = local_illumination(final_sdf, final_result.pos, dir, obj_col, inv_light_dir, light_col, ambient, intensity, spec_pow, shadow_pow, a_occ, final_result.steps, epsilon, max_dist, max_steps, dist);
        }
    }
    
    return final_col;
}