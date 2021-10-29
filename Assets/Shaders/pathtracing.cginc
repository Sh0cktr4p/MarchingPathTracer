#include "raymarching.cginc"

#define PI 3.1415926

float rand(float2 uv){
    return frac(sin(_Time.z + uv.x * uv.y + 12345 * uv.y / uv.x) * 43758.5453123);
}

float3 jitter(float3 dir, float phi, float sina, float cosa){
    float3 w = normalize(dir);
    float3 u = normalize(cross(w.yzx, w));
    float3 v = cross(w, u);
    return (u * cos(phi) + v * sin(phi)) * sina + w * cosa;
}


class PathTracer {
    Ray _view_ray;
    RayMarcher _ray_marcher;
    float3 _inv_light_dir;
    int _max_bounces;
    float2 uv;
    float3 _obj_col;
    float _intensity;
    float _roughness;
    float _focal_dist;
    float _aperture;
    float _lamp_rad;
    float3 _cam_fwd;
    float3 _cam_right;
    float3 _cam_up;
    int _random_seed;

    float get_random_value() {
        float random_value = rand(_uv * cam_up.xy + _Time.xz + _random_seed + 2.325876);
        _random_seed++;
        return random_value;
    }

    // Used to offset the starting ray to get a depth of field effect
    Ray simulate_lens() {
        float angle = get_random_value() * 2 * PI;
        float radius = sqrt(get_random_value());
        float2 offset = float2(cos(angle), sin(angle)) * radius * aperture;
        float aperture_area = PI * aperture * aperture;

        float4 focal_plane = float4(-_cam_fwd, dot(_cam_fwd, (_view_ray._pos + _cam_fwd * _focal_dist)));
        float dist_to_focal_plane = -(dot(_view_ray._pos, focal_plane.xyz) + focal_plane.w) / dot(_view_ray._dir, focal_plane.xyz);

        float3 focus_pos = _view_ray._pos + dist_to_focal_plane * _view_ray._dir;

        Ray distorted_ray = to_ray (
            _view_ray._pos + _cam_right * offset.x + _cam_up * offset.y,
            normalize(focus_pos - _view_ray._pos),
        );

        return distorted_ray;
    }

    float3 get_next_event_lighting(SDF light_sdf, LightingModel lighting_model, Ray in_ray, float3 pos, float3 nor) {
        // Get the ray to the closest point on the light source
        Ray next_event_ray = to_ray(pos, _ray_marcher.calc_normal(light_sdf, pos) * -1.);

        MarchResult next_event_march_result = _ray_marcher.march(sdf_with_lamp, next_event_ray);

        bool hit_light_source = _ray_marcher.hit_sdf(light_sdf, next_event_march_result);

        return _intensity * lighting_model.eval(n, in_ray._dir, next_event_ray._dir) * hit_light_source;
    }

    Ray reflect_ray(Ray in_ray, float3 pos) {

    }

    fixed4 trace_path(SDF scene_sdf, SDF light_sdf, LightingModel lighting_model) {
        Ray out_ray = simulate_lens();

        float3 color_mask = 1;
        float3 acc_color = 0;

        MarchResult result = _ray_marcher.march(scene_sdf, out_ray);

        if (!_ray_marcher.hit_something(result)) {
            return saturate(fixed4(acc_col, 1));
        }

        SDF combined_sdf = unite(scene_sdf, light_sdf);

        for (int i = 0; i < _max_bounces; i++) {
            Ray in_ray = out_ray;
            float3 nor = _ray_marcher.calc_normal(combined_sdf, result._pos);

            // Check whether we hit a light source, otherwise perform next event estimation
            if (_ray_marcher.hit_sdf(light_sdf, result)) {
                acc_col += col_mask * _intensity;
            }
            else {
                acc_col += col_mask * get_next_event_lighting(light_sdf, lighting_model, in_ray, result._pos, nor);
            }

            // Reflect ray
            out_ray = reflect_ray(in_ray, result._pos);
            col_mask *= lighting_model.eval(nor, in_ray, out_ray);

            // March on
            result = _ray_marcher.march(combined_sdf, out_ray);

            if (!_ray_marcher.hit_something(result)) {
                return saturate(fixed4(acc_col, 1));
            }
        }

        return saturate(fixed4(acc_col, 1));
    }
};

PathTracer getPathTracer(Ray view_ray, float3 inv_light_dir, float epsilon, float max_dist, int max_steps, int max_bounces, float2 uv, float3 obj_col, float intensity, float roughness,
        float focal_dist, float aperture, float3 cam_fwd, float3 cam_right, float3 cam_up) {
    PathTracer p = {
        view_ray,
        inv_light_dir,
        epsilon,
        max_dist,
        max_steps,
        max_bounces,
        uv,
        obj_col,
        intensity,
        roughness,
        focal_dist,
        aperture,
        cam_fwd,
        cam_right,
        cam_up,
        0,
    };
    return p;
}


fixed4 trace_path(SDF sdf, float3 pos, float3 dir, float3 inv_light_dir, float epsilon, float max_dist, int max_steps, int max_bounces, float2 uv, float3 obj_col, float intensity,
        float focal_dist, float aperture, float3 cam_fwd, float3 cam_right, float3 cam_up, float roughness){
    
    float3 cam_pos = pos;
    float3 curr_pos = pos;
    float3 curr_dir = dir;

    // Specify aperture offset
    float angle = rand(uv * cam_right.xz + _Time.xz + 2.325876) * 2.0 * PI;
    float radius = sqrt(rand(uv * cam_up.xy + _Time.zx * 4 + PI * PI));
    float2 offset = float2(cos(angle), sin(angle)) * radius * aperture;
    float aperture_area = PI * aperture * aperture;

    // Intensity has to be modified according to aperture area
    //intensity /= aperture_area;

    // Simulate camera lens
    float4 focal_plane = float4(-cam_fwd, dot(cam_fwd, (cam_pos + cam_fwd * focal_dist)));
    float dist_to_focal_plane = -(dot(cam_pos, focal_plane.xyz) + focal_plane.w) / dot(curr_dir, focal_plane.xyz);

    float3 focus_pos = cam_pos + dist_to_focal_plane * curr_dir;

    // Apply aperture offset
    curr_pos = curr_pos + cam_right * offset.x + cam_up * offset.y;
    curr_dir = normalize(focus_pos - curr_pos);


    float lamp_rad = 0.5;

    float3 col_mask = 1;
    float3 acc_col = 0;

    MarchResult result = march(sdf, curr_pos, curr_dir, epsilon, max_dist, max_steps);

    float dist = length(result.pos - curr_pos);

    if(dist < epsilon){
        return 0;//fixed4(1.0, 0.0, 1.0, 1.0);
    }

    SDF sdf_new = unite(sdf, translate(sphere(lamp_rad), cam_pos));
    
    for(int i = 0; i < max_bounces; i++){
        curr_pos = result.pos;
        
        float3 n = calcNormal(sdf_new, curr_pos, epsilon);
        
        float3 l0 = cam_pos - curr_pos;
		float cos_a_max = sqrt(1. - clamp(lamp_rad * lamp_rad / dot(l0, l0), 0., 1.));
		float cosa = lerp(cos_a_max, 1., rand(uv + _Time.xz + pos.xy));
		float3 l = jitter(l0, 2.*PI*rand(uv + _Time.xz + pos.zx), sqrt(1. - cosa*cosa), cosa);
        

        // Next-event estimation
        float ne_dist = april(sdf_new, curr_pos + l * epsilon * 2, l, epsilon, max_dist, max_steps) + 2 * epsilon;
        if(length(cam_pos - (curr_pos + l * ne_dist)) < lamp_rad + epsilon){
            float omega = 2.0 * PI * (1.0 - cos_a_max);
            //acc_col += (col_mask * obj_col * saturate(pow(dot(-reflect(curr_dir, n), l), 4)) * intensity * omega) / PI;
            //acc_col += (col_mask * obj_col * saturate(pow(dot(normalize(l - curr_dir), n), 4)) * intensity * omega) / PI;
            //acc_col += (col_mask * obj_col * intensity * omega * torrance_sparrow(n, l, curr_dir, 0.2)) / PI;
            acc_col += col_mask * intensity * omega * torrance_sparrow(n, l, curr_dir, roughness);
        }
        else if(i == 0) {
        }
            return length(cam_pos - (curr_pos + l * ne_dist));

        float3 nl = n * sign(-dot(n, curr_dir));

        float r2 = rand(uv + _Time.xy);

        col_mask *= torrance_sparrow(n, l, curr_dir, roughness);
        curr_dir = rand(uv + 3.3 * _Time.xy) < roughness ? jitter(nl, 2 * PI * rand(uv + _Time.xy + 1), sqrt(r2), sqrt(1 - r2)) : curr_dir;
        curr_dir = reflect(curr_dir, n);

        // Random light source hits
        if(length(cam_pos - curr_pos) < lamp_rad + epsilon){
            acc_col += col_mask * intensity;
            curr_dir = jitter(nl, 2 * PI * rand(uv + _Time.xy + 1), sqrt(r2), sqrt(1 - r2));
        }

        //col_mask *= obj_col * dot(l, reflect(curr_dir, n));

        result = march(sdf_new, curr_pos + curr_dir * epsilon * 2, curr_dir, epsilon, max_dist, max_steps);

        dist = length(result.pos - curr_pos);// + 2 * epsilon;

        if(dist < epsilon){
            return saturate(fixed4(acc_col, 1));
        }
    }

    return saturate(fixed4(acc_col, 1));
}