Shader "Custom/InfiniteSpheres"
{
    Properties
    {
        _ObjCol ("Object Color", Color) = (0.0, 1.0, 0.0, 1.0)
        _BackCol ("Background Color", Color) = (0.0, 0.0, 0.0, 1.0)
        _Roughness ("Roughness", Float) = 0.2
        _Epsilon ("Epsilon", Float) = 0.001
        _MaxDist ("Maximum Distance", Float) = 100
        _MaxSteps ("Maximum Steps", Int) = 100

        _A ("Parameter A", Float) = 0.5
        _B ("Parameter B", Float) = 10
        _C ("Parameter C", Float) = 10
    }
    SubShader
    {
        // No culling or depth
        Cull Off ZWrite Off ZTest Always

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"
            #include "UnityLightingCommon.cginc"
            #include "./pathtracing.cginc"

            struct appdata
            {
                float4 pos : POSITION;
                float2 tex : TEXCOORD0;
            };

            struct v2f
            {
                float4 pos : SV_POSITION;
                float2 tex : TEXCOORD0;
                float3 ray : TEXCOORD1;
            };

            sampler2D _MainTex;
            sampler2D _BackgroundTex;
            float4 _ObjCol;
            float4 _BackCol;
            float _Intensity;
            float _Roughness;

            float _Epsilon;
            float _MaxDist;
            float _MaxSteps;

            float4x4 _FrustumCornersMatrix;
            float4x4 _InvViewMatrix;
            float4 _MainTex_TexelSize;
            float3 _CameraPos;
            float3 _CameraFwd;
            float3 _CameraUp;
            float3 _CameraRight;
            int _FeedbackCounter;

            float _FocalDist;
            float _Aperture;

            float _A;
            float _B;
            float _C;

            v2f vert (appdata v)
            {
                v2f o;
                half index = v.pos.z;
                v.pos.z = 0.1;
                o.pos = UnityObjectToClipPos(v.pos);
                o.tex = v.tex;

                #if UNITY_UV_STARTS_AT_TOP
                if(_MainTex_TexelSize.y < 0){
                    o.tex.y = 1 - o.tex.y;
                }
                #endif

                o.ray = _FrustumCornersMatrix[(int)index].xyz;
                
                o.ray = mul(_InvViewMatrix, o.ray);

                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                fixed4 rec = tex2D(_MainTex, i.tex);
                fixed4 back = tex2D(_BackgroundTex, i.tex);
                back *= 0.0;
            

                //SDF scene_sdf = torus(float2(100, 30) * 2);
                SDF scene_sdf = menger_sponge(4, 320);  
                //SDF scene_sdf = menger_lod(4, 320, 8, _CameraPos, normalize(i.ray), _Epsilon, _MaxDist, _MaxSteps);
                
                scene_sdf = sphere_fold(scene_sdf, _A, _B);
                scene_sdf = repeat(scene_sdf, _C);
                
                class Custom: SDF {
                    float eval(float3 pos){
                        //pos = mirror_warp(pos, )
                        return sphere(1).eval(pos);
                    }
                };

                //sdf = sphere(1);

                //Custom scene_sdf;


                SDF light_sdf = sphere(20);
                //SDF light_sdf = rotateX(translate(torus(float2(200, 60)), float3(200, 0, 0)), 90);
                //light_sdf = translate(light_sdf, _CameraPos + _CameraUp * 30);
                //light_sdf = repeat(light_sdf, 200);

                //LightingModel lighting_model = torrance_sparrow(_Roughness, COPPER);
                LightingModel lighting_model = lambert(GOLD);


                int max_bounces = 3;
                float3 view_dir = normalize(i.ray);
                Ray view_ray = to_ray(_CameraPos, view_dir);

                PathTracer path_tracer = new_path_tracer(view_ray, _WorldSpaceLightPos0.xyz, _Epsilon, _MaxDist, _MaxSteps, max_bounces, i.tex, _Intensity * 40, _FocalDist * 0.00001, _Aperture, _CameraRight, _CameraUp, _CameraFwd);

                fixed4 res = saturate(path_tracer.trace_path(scene_sdf, light_sdf, lighting_model));
                //return fixed4(res.xyz * res.w + col.xyz * (1 - res.w), 1);
                

                return fixed4(1.0 * ((res.xyz * res.a + back.xyz * (1 - res.a)) + rec.xyz * _FeedbackCounter) / (_FeedbackCounter + 1), 1);
                //return fixed4((res.xyz * res.a + back.xyz * (1 - res.a) + rec.xyz * (_FeedbackCounter + 1)) / (_FeedbackCounter + 2), 1);

            }
            ENDCG
        }
    }
}
