#include "./quaternion.cginc"

interface SDF{
    float eval(float3 pos);
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++ PRIMITIVES +++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plane: SDF {
    float3 _nor;
    float _off;

    float eval(float3 pos){
        return dot(pos, _nor) + _off;
    }
};

Plane plane(float3 nor, float off){
    Plane p = {normalize(nor), off};
    return p;
}


class Sphere: SDF {
    float _rad;

    float eval(float3 pos){
        return length(pos) - _rad;
    }
};

Sphere sphere(float rad){
    Sphere s = {rad};
    return s;
}


class Ellipsoid: SDF {
    float3 _rad;

    float eval(float3 pos){
        float k0 = length(pos / _rad);
        float k1 = length(pos / (_rad * _rad));
        return k0*(k0-1.0)/k1;
    }
};

Ellipsoid ellipsoid(float3 rad){
    Ellipsoid e = {rad};
    return e;
}


class Box: SDF {
    float3 _dim;

    float eval(float3 pos){
        float3 d = abs(pos) - _dim;
        return length(max(d, 0.0)) + min(max(d.x, max(d.y, d.z)), 0.0);
    }
};

Box box(float3 dim){
    Box b = {dim};
    return b;
}


class Cylinder: SDF {
    float3 _dim;

    float eval(float3 pos){
        return length(pos.xz - _dim.xy) - _dim.z;
    }
};

Cylinder cylinder(float3 dim){
    Cylinder c = {dim};
    return c;
}


class Torus: SDF {
    float2 _thc;

    float eval(float3 pos){
        float2 q = float2(length(pos.xz) - _thc.x, pos.y);
        return length(q) - _thc.y;
    }
};

Torus torus(float2 thc){
    Torus t = {thc};
    return t;
}


class Tetrahedron: SDF {
    float _rad;

    float eval(float3 pos){
        float md = max(max(-pos.x - pos.y - pos.z, pos.x + pos.y - pos.z), max(-pos.x + pos.y + pos.z, pos.x - pos.y + pos.z));
        return (md - _rad) / 1.732050807;
    }
};

Tetrahedron tetrahedron(float rad){
    Tetrahedron t = {rad};
    return t;
}


class Octahedron: SDF {
    float _rad;

    float eval(float3 pos){
        pos = abs(pos);
        return (pos.x + pos.y + pos.z - _rad) * 0.57735027;
    }
};

Octahedron octahedron(float rad){
    Octahedron o = {rad};
    return o;
}


class Pyramid: SDF {
    float _dim;

    float eval(float3 pos){
        float m2 = _dim * _dim + 0.25;
    
        pos.xz = abs(pos.xz);
        pos.xz = (pos.z > pos.x) ? pos.zx : pos.xz;
        pos.xz -= 0.5;

        float3 q = float3( pos.z, _dim * pos.y - 0.5 * pos.x, _dim * pos.x + 0.5 * pos.y);
   
        float s = max(-q.x, 0.0);
        float t = clamp((q.y - 0.5 * pos.z) / (m2 + 0.25), 0.0, 1.0);
    
        float a = m2 * (q.x + s) * (q.x + s) + q.y * q.y;
        float b = m2 * (q.x + 0.5 * t) * (q.x + 0.5 * t) + (q.y - m2 * t) * (q.y - m2 * t);
    
        float d2 = min(q.y, -q.x * m2-q.y * 0.5) > 0.0 ? 0.0 : min(a, b);
    
        return sqrt((d2 + q.z * q.z) / m2) * sign(max(q.z, -pos.y));
    }
};

Pyramid pyramid(float dim){
    Pyramid p = {dim};
    return p;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++ COMBINATORS ++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SDF intersect(SDF fna, SDF fnb){
    class Intersect: SDF {
        float eval(float3 pos){
            return max(fna.eval(pos), fnb.eval(pos));
        }
    };

    Intersect i;

    return i;
}

SDF intersect_smooth(SDF fna, SDF fnb, float amt){
    class IntersectSmooth: SDF {
        float eval(float3 pos) {
            float da = fna.eval(pos);
            float db = fnb.eval(pos);
            float h = clamp(0.5 - 0.5 * (db - da) / amt, 0.0, 1.0);
            return lerp(db, da, h) + amt * h * (1.0 - h);
        }
    };

    IntersectSmooth i;
    return i;
}

SDF unite(SDF fna, SDF fnb){
    class Unite: SDF {
        float eval(float3 pos) {
            return min(fna.eval(pos), fnb.eval(pos));
        }
    };

    Unite u;
    return u;
}

SDF unite_smooth(SDF fna, SDF fnb, float amt){
    class UniteSmooth: SDF {
        float eval(float3 pos) {
            float da = fna.eval(pos);
            float db = fnb.eval(pos);
            float h = clamp(0.5 + 0.5 * (db - da) / amt, 0.0, 1.0);
            return lerp(db, da, h) - amt * h * (1.0-h);
        }
    };

    UniteSmooth u;
    return u;
}

SDF subtract(SDF fna, SDF fnb){
    class Subtract: SDF {
        float eval(float3 pos) {
            return max(fna.eval(pos), -fnb.eval(pos));
        }
    };

    Subtract s;
    return s;
}

SDF subtract_smooth(SDF fna, SDF fnb, float amt){
    class SubtractSmooth: SDF {
        float eval(float3 pos) {
            float da = fna.eval(pos);
            float db = fnb.eval(pos);
            float h = clamp( 0.5 - 0.5 * (db + da) / amt, 0.0, 1.0);
            return lerp(db, -da, h) + amt * h * (1.0 - h);
        }
    };

    SubtractSmooth s;
    return s;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++ TRANSFORMATORS +++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
float3 translate_warp(float3 pos, float3 amt){
    return pos - amt;
}

SDF translate(SDF sdf, float3 amt){
    class Translate: SDF {
        float eval(float3 pos){
            return sdf.eval(translate_warp(pos, amt));
        }
    };

    Translate t;
    return t;
}


float3 rotate_warp(float3 pos, quaternion rot){
    return qrot(rot, pos);
}

SDF rotate(SDF sdf, quaternion rot){
    class Rotate: SDF {
        float eval(float3 pos){
            return sdf.eval(rotate_warp(pos, rot));
        }
    };

    Rotate r;

    return r;
}


float3 rotateX_warp(float3 pos, float deg){
    float s, c;
    sincos(deg * 3.1415926 * -1 / 180, s, c);
    return float3(pos.x, c * pos.y - s * pos.z, s * pos.y + c * pos.z);
}

SDF rotateX(SDF sdf, float deg){
    class RotateX: SDF {
        float eval(float3 pos){
            return sdf.eval(rotateX_warp(pos, deg));
        }
    };

    RotateX r;
    return r;
}


float3 rotateY_warp(float3 pos, float deg){
    float s, c;
    sincos(deg * 3.1415926 * -1 / 180, s, c);
    return float3(c * pos.x + s * pos.z, pos.y, c * pos.z - s * pos.x);
}

SDF rotateY(SDF sdf, float deg){
    class RotateY: SDF {
        float eval(float3 pos){
            return sdf.eval(rotateY_warp(pos, deg));
        }
    };

    RotateY r;
    return r;
}


float3 rotateZ_warp(float3 pos, float deg){
    float s, c;
    sincos(deg * 3.1415926 * -1 / 180, s, c);
    return float3(c * pos.x - s * pos.y, s * pos.x + c * pos.y, pos.z);
}

SDF rotateZ(SDF sdf, float deg){
    class RotateZ: SDF {
        float eval(float3 pos){
            return sdf.eval(rotateZ_warp(pos, deg));
        }
    };

    RotateZ r;
    return r;
}


SDF scale(SDF sdf, float amt){
    class Scale: SDF {
        float eval(float3 pos) {
            return sdf.eval(pos / amt) * amt;
        }
    };

    Scale s;
    return s;
}


float3 elongate_warp(float3 pos, float3 amt){
    float3 q = abs(pos) - amt;
    return max(q, 0.0) + min(max(q.x, max(q.y, q.z)), 0.0);
}

SDF elongate(SDF sdf, float3 amt){
    class Elongate: SDF {
        float eval(float3 pos) {
            float3 q = abs(pos) - amt;
            return sdf.eval(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
        }
    };

    Elongate e;
    return e;
}


float3 repeat_warp(float3 pos, float gap){
    return (frac((pos + 0.5 * gap) / gap) - 0.5) * gap;
}

SDF repeat(SDF sdf, float gap){
    class Repeat: SDF {
        float eval(float3 pos){
            return sdf.eval(repeat_warp(pos, gap));
        }
    };

    Repeat r;
    return r;
}


float3 repeat_finite_warp(float3 pos, float gap, float3 amt){
    return pos - gap * clamp(round(pos / gap), -amt, amt);
}

SDF repeat_finite(SDF sdf, float gap, float3 amt){
    class RepeatFinite: SDF {
        float eval(float3 pos) {
            return sdf.eval(repeat_finite_warp(pos, gap, amt));
        }
    };
}


float3 mirror_warp(float3 pos, float3 nor, float off){
    return pos - 2.0 * min(0.0, dot(pos, nor) - off) * nor;
}

SDF mirror(SDF sdf, float3 nor, float off){
    class Mirror: SDF {
        float eval(float3 pos){
            return sdf.eval(mirror_warp(pos, nor, off));
        }
    };

    Mirror m;
    return m;
}

float reflectXY(float3 pos, SDF sdf, float off){
    return sdf.eval(float3(pos.x, pos.y, abs(pos.z + off) - off));
}

float reflectXZ(float3 pos, SDF sdf, float off){
    return sdf.eval(float3(pos.x, abs(pos.y + off) - off, pos.z));
}

float reflectYZ(float3 pos, SDF sdf, float off){
    return sdf.eval(float3(abs(pos.x + off) - off, pos.y, pos.z));
}

float round(float3 pos, SDF sdf, float rad){
    return sdf.eval(pos) - rad;
}

float onion(float3 pos, SDF sdf, float thc){
    return abs(sdf.eval(pos)) - thc;
}


SDF sin_displace(SDF sdf, float frq, float amt){
    class SinDisplace: SDF {
        float eval(float3 pos){
            float d1 = sdf.eval(pos);
            float d2 = sin(frq * pos.x) * sin(frq * pos.y) * sin(frq * pos.z) * amt;
            return d1 + d2;
        }
    };

    SinDisplace s;
    return s;
}


float3 twistY_warp(float3 pos, float amt){
    float s, c;
    sincos(pos.y * amt * 3.1415926 / 180, s, c);
    float2x2 m = float2x2(c, -s, s, c);
    return float3(mul(m, pos.xz), pos.y);
}

SDF twistY(SDF sdf, float amt){
    class TwistY: SDF {
        float eval(float3 pos){
            return sdf.eval(twistY_warp(pos, amt));
        }
    };

    TwistY t;
    return t;
}


float3 sphere_fold_warp(float3 pos, float rMin, float rMax){
    float r2 = dot(pos, pos);
    return pos * max(rMax / max(rMin, r2), 1.0);
}

SDF sphere_fold(SDF sdf, float mr2, float rd2){
    class SphereFold: SDF {
        float eval(float3 pos){
            float r2 = dot(pos, pos);
            float factor = max(rd2 / max(mr2, r2), 1.0);
            //float factor = min(mr2 / min(rd2, r2), 1.0);
            return sdf.eval(pos * factor) / factor;
            //return sdf.eval(sphere_fold_warp(pos, mr2, rd2));
        }
    };

    SphereFold s;
    return s;
}


float3 box_fold_warp(float3 pos, float3 lim) {
	return clamp(pos, -lim, lim) * 2.0 - pos;
}

SDF box_fold(SDF sdf, float3 lim){
    class BoxFold: SDF {
        float eval(float3 pos){
            return sdf.eval(box_fold_warp(pos, lim));
        }
    };

    BoxFold b;
    return b;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++ COMPLEX FUNCTIONS ++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SDF sierpinsky_tetrahedron(int lev, float scl){
    class Sierpinsky: SDF{
        float eval(float3 pos){
            float mirror_dist = -sqrt(2) / 2;
            float3 mirror_normal_a = normalize(float3(1, 1, 0));
            float3 mirror_normal_b = mirror_normal_a.xzy;
            float3 mirror_normal_c = mirror_normal_a.zxy;

            float scale = scl;

            for(int i = 0; i < lev; i++){
                pos = translate_warp(pos, scale * 0.5);
                pos = mirror_warp(pos, mirror_normal_a, mirror_dist * scale);
                pos = mirror_warp(pos, mirror_normal_b, mirror_dist * scale);
                pos = mirror_warp(pos, mirror_normal_c, mirror_dist * scale);

                scale *= 0.5;
            }

            return tetrahedron(scale).eval(pos);
        }
    };

    Sierpinsky s;
    return s;
}

float billo_march(SDF sdf, float3 pos, float3 dir, float epsilon, float max_dist, int max_steps){
    float dist = epsilon;
    float min_dist = 1.0;

    float3 curr_pos = pos + dir * dist;

    for(int i = 0; i < max_steps; i++){
        float new_dist = sdf.eval(curr_pos);
        dist += new_dist;
        curr_pos = pos + dir * dist;

        if(new_dist < epsilon){
            return dist;
        }

        if(dist > max_dist){
            return dist;
        }
    }

    return dist;
}


SDF menger_sponge(int lev, float scl){
    class Menger: SDF {
        float eval(float3 pos){
            float3 mirror_normal_a = normalize(float3(-1, 0, 1));
            float3 mirror_normal_b = normalize(float3(-1, 1, 0));
            float3 mirror_normal_c = normalize(float3(1, 0, 0));
            float scale = scl;

            float pulse = (max(sin(_Time.w*5), 0.5) - 0.5) * 2 * 0;

            for(int i = 0; i < lev; i++){
                scale /= 3;

                pos = translate_warp(pos, scale * 2);
                pos = translate_warp(pos, scale * float3(1, 0, 0) * 0.2);
                pos = abs(pos + scale * 2) - scale * 2;
                pos = mirror_warp(pos, rotateZ_warp(mirror_normal_a, 5+_Time.w*0), 0);
                pos = mirror_warp(pos, rotateX_warp(mirror_normal_b, 5+_Time.w*0), 0);
                pos = mirror_warp(pos, rotateY_warp(mirror_normal_c, 5+_Time.w*0), -scale);
            }

            return box(scale).eval(pos) - 0.25 * scale;//- pulse * 0.8 * scale;
        }
    };

    Menger m;
    return m;
}


SDF menger_sponge_classic(int lev, float scl) {
    class Menger : SDF {
        float eval(float3 pos) {
            float3 mirror_normal_a = normalize(float3(-1, 0, 1));
            float3 mirror_normal_b = normalize(float3(-1, 1, 0));
            float3 mirror_normal_c = normalize(float3(1, 0, 0));
            float scale = scl;

            float pulse = (max(sin(_Time.w * 5), 0.5) - 0.5) * 2 * 0;

            for (int i = 0; i < lev; i++) {
                scale /= 3;

                pos = translate_warp(pos, scale * 2);
                //pos = translate_warp(pos, scale * float3(1, 0, 0) * 0.2);
                pos = abs(pos + scale * 2) - scale * 2;
                pos = mirror_warp(pos, mirror_normal_a, 0);
                pos = mirror_warp(pos, mirror_normal_b, 0);
                pos = mirror_warp(pos, mirror_normal_c, -scale);
            }

            return box(scale).eval(pos) -0.25 * scale;//- pulse * 0.8 * scale;
        }
    };

    Menger m;
    return m;
}

SDF menger_lod(int lev, float scl, float lim, float3 pps, float3 dir, float eps, float mxd, int mxs){
    class MengerLOD: SDF {
        float eval(float3 pos) {
            float dist = billo_march(menger_sponge(lev, scl), pps, dir, eps, mxd, mxs);

            if(dist + 1 > mxd){
                return mxd;
            }

            int level = min(lev, log10(scl / (lim * dist)) / log10(3));
            
            return billo_march(menger_sponge(level, scl), pos, dir, eps, mxd, mxs);
            //return menger_sponge(level, scl).eval(pos);
        }
    };

    MengerLOD m;
    return m;
}

SDF mandelbox(int lev, float scl, float mr2, float rd2, float lim){
    class Mandelbox: SDF {
        float eval(float3 pos){
	        float3 offset = pos;
	        float dr = 1.0;
            float scale = scl;

	        for (int n = 0; n < lev; n++) {
		        pos = box_fold_warp(pos, lim);            // Reflect
		        pos = sphere_fold_warp(pos, mr2, rd2);    // Sphere Inversion
                dr = sphere_fold_warp(dr, mr2, rd2).x;
 		
                pos = scale * pos + offset;  // Scale & Translate
                dr = dr * abs(scale) + 1.0;
	        }

	        float r = length(pos);
	        return r / abs(dr);
        }
    };

    Mandelbox m;
    return m;
}

SDF warped_menger_sponge(int lev, int wct, float scl){
    class WarpedMengerSponge: SDF {
        float eval(float3 pos){
            float factor = 1;

            float rMax = scl * 64;
            float rMin = scl * 30;
            float scf = 1;

            for(int i = 0; i < wct; i++){
                pos = box_fold_warp(pos, scl * scf);
                scf *= 1.;
            }
                float r2 = dot(pos, pos);
                factor *= max(rMax / max(rMin, r2), 1.0);
                rMax /= 4;
                rMin /= 4;
                pos *= factor;
            

            return menger_sponge(lev, scl).eval(pos) / factor;
        }
    };

    WarpedMengerSponge w;
    return w;
}



float4 sphere_test(float4 z)
{
   float r2 = dot(z.xyz,z.xyz);
   if(r2<2.0)
      z*=(1.0/r2);
   else   z*=0.5;

   return z;

}
float3 box_test(float3 z){
   return clamp(z, -1.0, 1.0) * 2.0 - z;
}

float DE0(float3 pos, float3 from){
   float3 z=pos-from;
   float r=dot(pos-from,pos-from)*pow(length(z),2.0);
   return (1.0-smoothstep(0.0,0.01,r))*0.01;
}

float DE2(float3 pos, float3 params){
    float4 scale = -20*0.272321;
    float4 p = float4(pos,1.0), p0 = p;  
    float4 c = float4(params,0.5)-0.5; // param = 0..1

    for (float i=0;i<10; i++)
    {
        p.xyz = box_test(p.xyz);
        p = sphere_test(p);
        p = p * scale + c;
    }

    return length(p.xyz)/p.w;
}


SDF trippy(float3 from, float3 params){
    class Trippy: SDF {
        float eval(float3 pos){
            float d0 = DE0(pos, from);   
            float d2 = DE2(pos, params);

            return max(d0, d2);
        }
    };

    Trippy t;
    return t;
}






/*
float sdf_menger(float3 pos){
    float scalingFactor = 81.0;
    int limit = 1;
    for(int i = 0; i < limit; i++){
        pos = sdf_translate(pos, scalingFactor * 2 / 3.0);
        scalingFactor /= 3;
        pos = abs(pos + scalingFactor * 2) - scalingFactor * 2;
        pos = sdf_mirror(pos, rotY(normalize(float3(-1, 0, 1)), _Time.z * 0), 0);
        pos = sdf_mirror(pos, rotX(rotZ(rotY(normalize(float3(-1, 1, 0)), _Time.z * 5), _Time.z * 4.5), _Time.z * 2.1), 0);
        pos = sdf_mirror(pos, rotY(float3(1, 0, 0), _Time.z * 0), -scalingFactor);
    }
    return sdf_box(pos, scalingFactor);
}

float3 sdf_foldBox(float3 pos, float3 r){
    return clamp(pos, -r, r) * 2 - pos;
}

float3 sdf_foldSphere(float3 pos, float rMin, float rMax){
    float r2 = dot(pos, pos);
    return pos * max(rMax / max(rMin, r2), 1.0);
}


float sdf_infCross(float3 p, float r) {
    float3 q = p * p;
    return sqrt(min(min(q.x + q.y, q.x + q.z), q.y + q.z)) - r;
}

float sdf_infiniteSpheres(float3 pos){
    return sdf_sphere(frac(pos / 5.0) * 5.0 - 2.5);
}



float sdf_mandelbox(float3 pos){
    float scalingFactor = 8;
    int limit = 4;
    for(int i = 0; i < limit; i++){
        //pos = sdf_foldBox(pos, _SinTime.z);
        pos = sdf_foldSphere(pos, _SinTime.y, _SinTime.z);
        scalingFactor /= 2;
        
    }
    return sdf_box(pos, 6);
}

float sdf_mandelbulb(float3 pos){
    float3 z = pos;
    float dr = 1;
    float r;
    float power = 2;
    for(int i = 0; i < 15; i++){
        r = length(z);
        if(r > 2) break;
        float theta = acos(z.z / r) * power;
        float phi = atan2(z.y, z.x) * power;
        float zr = pow(r, power);
        dr = pow(r, power - 1) * power * dr + 1;
        z = zr * float3(sin(theta) * cos(phi), sin(phi) * sin(theta), cos(theta));
        z += pos; 
    }
    return 0.5 * log(r) * r / dr;
}

float sdf_test_cubes(float3 pos){
    return sdf_union(sdf_union(sdf_union(sdf_box(pos, float3(1, 0.25, 1)), 
                                        sdf_box(pos + float3(0.5, -0.0, 0.75), float3(1, 0.75, 0.75))),
                                sdf_union(sdf_box(pos + float3(-0.5, 0, 0.25), float3(0.5, 0.5, 0.75)),
                                        sdf_box(pos + float3(0.375, -0.25, 0.25), float3(0.125, 0.25, 0.5)))),
                    sdf_union(sdf_union(sdf_box(pos + float3(0.5, -0.0, 0.5), float3(0.25, 1, 0.25)), 
                                        sdf_box(pos + float3(0.75, 0.125, -0.5), float3(1.25, 0.125, 1.25))),
                                sdf_union(sdf_box(pos + float3(1.75, 0, 0.5), float3(0.5, 0.5, 0.75)),
                                        sdf_box(pos + float3(2.125, -0.0, -0.25), float3(0.375, 0.25, 0.5)))));
}

float sdf_test(float3 pos){
    return sdf_union(sdf_test_cubes(sdf_mirror(
                        sdf_mirror(pos, normalize(-float3(1, 0, 1)), -length(float3(0.5, 0, 0.5))),
                        normalize(float3(1, 0, -1)), -length(float3(1.25, 0, 1.25)))),
                        sdf_box(pos - float3(-0.75, 0, 1.75), float3(0.25, 1, 0.25)));
}*/