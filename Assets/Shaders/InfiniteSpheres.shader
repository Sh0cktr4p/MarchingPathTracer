Shader "Custom/InfiniteSpheres"
{
    Properties
    {
        _ObjCol ("Object Color", Color) = (0.0, 1.0, 0.0, 1.0)
        _BackCol ("Background Color", Color) = (0.0, 0.0, 0.0, 1.0)
        _Ambient ("Ambient", Range(0, 1)) = 0.1
        _SpecPow ("Specular Power", Float) = 1.
        _ShadowPow ("Shadow Power", Int) = 1
        _AmbOcc ("Ambient Occlusion", Float) = 200
        _Epsilon ("Epsilon", Float) = 0.001
        _MaxDist ("Maximum Distance", Float) = 100
        _MaxSteps ("Maximum Steps", Int) = 100
        //_FocalDist ("Focal Distance", Float) = 10
        //_Aperture ("Aperture", Float) = 1

        _MandelboxIts("Mandelbox Iterations", Int) = 10
        _MandelboxScl("Mandelbox Scale", Float) = 1
        _MandelboxMnr("Mandelbox Min Radius", Float) = 1
        _MandelboxMxr("Mandelbox Max Radius", Float) = 1
        _MandelboxLim("Mandelbox Box Limit", Float) = 1
        _ParamVector("Parameter Vector", Vector) = (0, 0, 0)
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
            //#include "./quaternion.cginc"
            //#include "./sdf.cginc"
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
            float4 _ObjCol;
            float4 _BackCol;
            float _Ambient;
            float _Intensity;
            float _SpecPow;
            int _ShadowPow;
            float _AmbOcc;
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

            int _MandelboxIts;
            float _MandelboxScl;
            float _MandelboxMnr;
            float _MandelboxMxr;
            float _MandelboxLim;
            float3 _ParamVector;

            float _FocalDist;
            float _Aperture;

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
                fixed4 col = tex2D(_MainTex, i.tex);
            
                SDF sdf = menger_sponge(4, 320);  
                //sdf = sphere_fold(sdf, _MandelboxMnr, _MandelboxMxr);
                //sdf = repeat(sdf, _MandelboxScl);
                
                /*class Custom: SDF {
                    float eval(float3 pos){
                        //pos = mirror_warp(pos, )
                        return sphere(1).eval(pos);
                    }
                };*/

                //sdf = sphere(1);

                //Custom sdf;

                //fixed4 res = menger_lod_march(_CameraPos, normalize(i.ray), _WorldSpaceLightPos0.xyz, _BackCol, float3(0.0, 0.2, 1.0), _LightColor0.xyz, _Ambient, _Intensity, _SpecPow, _ShadowPow, _AmbOcc, _Epsilon, _MaxDist, _MaxSteps, 12, 320, _MandelboxMnr, _MandelboxMxr, _MandelboxScl, 0.003);
            
                //fixed4 res = fixed4(0.0, 1.0, 1.0, 1.0);

                //fixed4 res = blinn_phong_raymarch(sdf, _CameraPos, normalize(i.ray), _WorldSpaceLightPos0.xyz, _BackCol, _ObjCol, _LightColor0.xyz, _Ambient, _Intensity, _SpecPow, _ShadowPow, _AmbOcc, _Epsilon, _MaxDist, _MaxSteps);
                fixed4 res = trace_path(sdf, _CameraPos, normalize(i.ray), _WorldSpaceLightPos0.xyz, _Epsilon, _MaxDist, _MaxSteps, 3, i.tex, float3(0.8, 0.0, 1.0), _Intensity * 40, _FocalDist, _Aperture, _CameraFwd, _CameraRight, _CameraUp, _SpecPow);

                
                //return fixed4(res.xyz * res.w + col.xyz * (1 - res.w), 1);
                return fixed4((res.xyz * res.w + col.xyz * col.w * _FeedbackCounter) / (_FeedbackCounter + 1), 1);
            }
            ENDCG
        }
    }
}
