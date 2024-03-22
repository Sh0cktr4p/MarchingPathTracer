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

// return a vector that is scattered in the hemisphere around nor. Angles in Radians
// nor: surface normal
// theta: azimuth
// phi: zenith
float3 scatter(float3 nor, float theta, float phi) {
    float3 scatter_vec = float3(sin(theta) * sin(phi), cos(phi), cos(theta) * sin(phi));

   // int over_surface = dot(nor, scatter_vec) > 0;
    //return scatter_vec;

    bool use_x = nor.x < 0.707 && nor.x > -0.707;
    float3 pseudo_bitangent = lerp(float3(0, 1, 0), float3(1, 0, 0), use_x);
    float3 tangent = normalize(cross(nor, pseudo_bitangent));
    float3 binormal = cross(nor, tangent);
    return scatter_vec.x * tangent + scatter_vec.y * nor + scatter_vec.z * binormal;
    //float3 tangent_space_normal = float3(0, 1, 0);
    //float3 rot_axis = normalize(cross(scatter_vec, nor));
    //float rot_axis = normalize(cross(nor, -float3(0, 1, 0)));



    //float angle = acos(dot(norm );
    //quaternion q = quat(rot_axis, angle);

    //return normalize(qrot(q, scatter_vec));
}

class PathTracer {
    Ray _view_ray;
    RayMarcher _ray_marcher;
    float3 _inv_light_dir;
    int _max_bounces;
    float2 _uv;
    float _intensity;
    float _focal_dist;
    float _aperture;
    float3 _cam_right;
    float3 _cam_up;
    float3 _cam_fwd;
    int _random_seed;

    float get_random_value() {
        float random_value = rand(_uv + sin(_Time.xz) + sin(_Time.y + 0.2) * 1.353 + _random_seed + 2.325876);
        _random_seed++;
        return random_value;
    }

    // Used to offset the starting ray to get a depth of field effect
    Ray simulate_lens() {
        float angle = get_random_value() * 2 * PI;
        float radius = sqrt(get_random_value());
        float2 offset = float2(cos(angle), sin(angle)) * radius * _aperture;
        float aperture_area = PI * _aperture * _aperture;

        float4 focal_plane = float4(-_cam_fwd, dot(_cam_fwd, (_view_ray._pos + _cam_fwd * _focal_dist)));
        float dist_to_focal_plane = -(dot(_view_ray._pos, focal_plane.xyz) + focal_plane.w) / dot(_view_ray._dir, focal_plane.xyz);

        float3 focus_pos = _view_ray._pos + dist_to_focal_plane * _view_ray._dir;

        Ray distorted_ray = to_ray (
            _view_ray._pos + _cam_right * offset.x + _cam_up * offset.y,
            normalize(focus_pos - _view_ray._pos)
        );

        return distorted_ray;
    }

    float3 get_next_event_lighting(SDF combined_sdf, SDF light_sdf, LightingModel lighting_model, Ray in_ray, float3 pos, float3 nor, float prev_dist) {
        // Get the ray to the closest point on the light source
        Ray next_event_ray = to_ray(pos, _ray_marcher.calc_normal(light_sdf, pos) * -1.);

        MarchResult next_event_march_result = _ray_marcher.march(combined_sdf, next_event_ray);

        bool hit_light_source = _ray_marcher.hit_sdf(light_sdf, next_event_march_result);
        float dist = prev_dist + length(next_event_march_result._pos - pos);

        return lighting_model.eval(nor, in_ray._dir, next_event_ray._dir) * hit_light_source * _intensity / (dist * dist);
    }

    Ray reflect_ray(Ray in_ray, float3 pos, float3 nor) {
        //return to_ray(pos, reflect(in_ray._dir, nor));
        float3 dir = scatter(nor, get_random_value() * 2 * PI, get_random_value() * 0.5 * PI);
        return to_ray(pos + 0.001 * dir, dir);
    }

    fixed4 trace_path(SDF scene_sdf, SDF light_sdf, LightingModel lighting_model) {
        Ray out_ray = simulate_lens();

        float3 color_mask = 1;
        float3 acc_color = 0;
        float acc_dist = 0;

        SDF combined_sdf = unite(scene_sdf, light_sdf);
        MarchResult result = _ray_marcher.march(combined_sdf, out_ray);

        if (!_ray_marcher.hit_something(result)) {
            return saturate(fixed4(acc_color, 0));
        }


        for (int i = 0; i < _max_bounces; i++) {
            Ray in_ray = out_ray;
            acc_dist += length(result._pos - out_ray._pos);
            float3 nor = _ray_marcher.calc_normal(combined_sdf, result._pos);
            float cos_theta = dot(nor, in_ray._dir);

            // Check whether we hit a light source, otherwise perform next event estimation
            if (_ray_marcher.hit_sdf(light_sdf, result)) {
                acc_color += color_mask * _intensity;// / (acc_dist * acc_dist);
            }
            else {
                acc_color += cos_theta * color_mask * get_next_event_lighting(combined_sdf, light_sdf, lighting_model, in_ray, result._pos, nor, acc_dist);
            }

            // Reflect ray
            out_ray = reflect_ray(in_ray, result._pos, nor);
            color_mask *= cos_theta * lighting_model.eval(nor, in_ray._dir, out_ray._dir);
            
            //if (length(out_ray._dir) < 0.5) return fixed4(1, 0, 0, 1);
            //return fixed4(dot(out_ray._dir, nor).xxx, 1);
            //return fixed4(nor - out_ray._dir, 1);
            // March on
            result = _ray_marcher.march(combined_sdf, out_ray);

            if (!_ray_marcher.hit_something(result)) {
                if (_ray_marcher.out_of_steps(result)) {
                    return fixed4(i == 0, i == 1, i == 2, 1) * fixed4(acc_color, 1);
                }
                else {
                    acc_color += color_mask * _intensity * 0.01;// *max(0.01, dot(out_ray._dir, float3(0, 1, 0)));
                }
                return saturate(fixed4(acc_color, 1));
            }
        }

        return saturate(fixed4(acc_color, 1));
    }
};

PathTracer new_path_tracer(Ray view_ray, float3 inv_light_dir, float epsilon, float max_dist, int max_steps, int max_bounces, float2 uv, float intensity,
        float focal_dist, float aperture, float3 cam_right, float3 cam_up, float3 cam_fwd) {
    RayMarcher ray_marcher = new_ray_marcher(epsilon, max_dist, max_steps);

    PathTracer p = {
        view_ray,
        ray_marcher,
        inv_light_dir,
        max_bounces,
        uv,
        intensity,
        focal_dist,
        aperture,
        cam_right,
        cam_up,
        cam_fwd,
        0
    };
    return p;
}
