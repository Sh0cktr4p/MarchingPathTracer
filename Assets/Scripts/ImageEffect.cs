using System.Linq;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
[RequireComponent(typeof(Camera))]
[AddComponentMenu("Effects/Raymarch (Generic)")]
public class ImageEffect : MonoBehaviour
{
    private RenderTexture tempTexture;

    [SerializeField]
    private Material effectMaterial;

    private Camera currentCamera;

    public Camera CurrentCamera
    {
        get
        {
            if(!currentCamera)
            {
                currentCamera = GetComponent<Camera>();
            }

            return currentCamera;
        }
    }

    public float intensity = 0f;
    public float aperture = 0.1f;
    public float focalDist = 30f;

    public bool feedback = false;

    private int feedbackCounter = 0;


    private Matrix4x4 GetFrustumCorners(Camera cam)
    {   
        float camFov = cam.fieldOfView;
        float camAspect = cam.aspect;

        Matrix4x4 frustumCorners = Matrix4x4.identity;

        float fovWHalf = camFov * 0.5f;

        float tan_fov = Mathf.Tan(fovWHalf * Mathf.Deg2Rad);

        Vector3 toRight = Vector3.right * tan_fov * camAspect;
        Vector3 toTop = Vector3.up * tan_fov;

        Vector3 topLeft = (-Vector3.forward - toRight + toTop).normalized;
        Vector3 topRight = (-Vector3.forward + toRight + toTop).normalized;
        Vector3 bottomRight = (-Vector3.forward + toRight - toTop).normalized;
        Vector3 bottomLeft = (-Vector3.forward - toRight - toTop).normalized;

        frustumCorners.SetRow(0, topLeft);
        frustumCorners.SetRow(1, topRight);
        frustumCorners.SetRow(2, bottomRight);
        frustumCorners.SetRow(3, bottomLeft);

        return frustumCorners;
    }

    private static void CustomBlit(RenderTexture src, RenderTexture dst, Material fxMat, int passNr){
        RenderTexture.active = dst;

        fxMat.SetTexture("_MainTex", src);

        GL.PushMatrix();
        GL.LoadOrtho();

        fxMat.SetPass(passNr);

        GL.Begin(GL.QUADS);

        GL.MultiTexCoord2(0, 0.0f, 0.0f);
        GL.Vertex3(0.0f, 0.0f, 3.0f);

        GL.MultiTexCoord2(0, 1.0f, 0.0f);
        GL.Vertex3(1.0f, 0.0f, 2.0f);

        GL.MultiTexCoord2(0, 1.0f, 1.0f);
        GL.Vertex3(1.0f, 1.0f, 1.0f);

        GL.MultiTexCoord2(0, 0.0f, 1.0f);
        GL.Vertex3(0.0f, 1.0f, 0.0f);

        GL.End();
        GL.PopMatrix();
    }

    [ImageEffectOpaque]
    void OnRenderImage(RenderTexture src, RenderTexture dst)
    {
        if(!effectMaterial)
        {
            Graphics.Blit(src, dst);
            return;
        }

        if(!tempTexture)
        {
            tempTexture = RenderTexture.GetTemporary(src.descriptor);
            //Graphics.Blit(src, tempTexture);
        }

        effectMaterial.SetMatrix("_FrustumCornersMatrix", GetFrustumCorners(CurrentCamera));
        effectMaterial.SetMatrix("_InvViewMatrix", CurrentCamera.cameraToWorldMatrix);
        effectMaterial.SetVector("_CameraPos", transform.position);
        effectMaterial.SetFloat("_MaxEpsilon", CurrentCamera.fieldOfView / CurrentCamera.pixelWidth);
        effectMaterial.SetFloat("_Intensity", intensity);
        effectMaterial.SetFloat("_FocalDist", focalDist);
        effectMaterial.SetFloat("_Aperture", aperture);
        effectMaterial.SetVector("_CameraFwd", CurrentCamera.transform.forward);
        effectMaterial.SetVector("_CameraRight", CurrentCamera.transform.right);
        effectMaterial.SetVector("_CameraUp", CurrentCamera.transform.up);
        effectMaterial.SetTexture("_BackgroundTex", src);

        Debug.Log(feedbackCounter);
        if(feedback)
        {
            effectMaterial.SetInt("_FeedbackCounter", feedbackCounter);
            CustomBlit(tempTexture, dst, effectMaterial, 0);
            Graphics.Blit(dst, tempTexture);
            feedbackCounter += 1;
        }
        else
        {
            Debug.Log("Resetting feedback counter");
            feedbackCounter = 0;
            effectMaterial.SetInt("_FeedbackCounter", feedbackCounter);
            CustomBlit(src, dst, effectMaterial, 0);
            Graphics.Blit(src, tempTexture);
        }

    }
}
