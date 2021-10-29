using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public struct LightProps{
    public Vector3 position;
    public Vector3 color;
    public float intensity;

    public LightProps(Vector3 pos, Vector3 col, float lin) {
        position = pos;
        color = col;
        intensity = lin;
    }
}

public class Light : MonoBehaviour
{
    public Color lightColor;
    public float lightIntensity;

    public LightProps Props{
        get {
            return new LightProps(transform.position, (Vector4)lightColor, lightIntensity);
        }
    }
}
