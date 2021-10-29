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
        return to_result(ray._pos, min_dist, _max_steps);
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
