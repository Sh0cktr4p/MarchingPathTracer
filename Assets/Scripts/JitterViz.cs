using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class JitterViz : MonoBehaviour
{
    public float phi;
    public float sina;
    public float cosa;
    public float theta;

    public Vector3 dir; 

    private Vector3 u;
    private Vector3 v;
    private Vector3 w;

    // Start is called before the first frame update
    void Start()
    {
        w = dir.normalized;
        u = Vector3.Cross(new Vector3(w.y, w.z, w.x), w).normalized;
        v = Vector3.Cross(w, u);
    }

    Vector3 Jitter(Vector3 d, float p, float s, float c)
    {
        Vector3 w = d.normalized;
        Vector3 u = Vector3.Cross(new Vector3(w.y, w.z, w.x), w).normalized;
        Vector3 v = Vector3.Cross(w, u);
        return u * Mathf.Cos(p) + v * Mathf.Sin(p) * s + w * c;
    }

    // Update is called once per frame
    void Update()
    {
        Vector3 newDir = u * Mathf.Cos(phi) + v * Mathf.Sin(phi) * sina + w * cosa;
        Vector3 newDir2 = u * Mathf.Cos(phi) + v * Mathf.Sin(phi) * Mathf.Sin(theta) + w * Mathf.Cos(theta);

        for(int i = 0; i < 100; i++)
        {
            float r1 = Random.value;
            float r2 = Random.value;
            float n1 = Mathf.Sqrt(-2 * Mathf.Log(r1)) * Mathf.Cos(2 * Mathf.PI * r2);
            float n2 = Mathf.Sqrt(-2 * Mathf.Log(r1)) * Mathf.Sin(2 * Mathf.PI * r2);
            Vector3 rayDir = Jitter(dir, 2 * Mathf.PI * n1, Mathf.Sqrt(n2), Mathf.Sqrt(1 - n2));
            Debug.DrawRay(Vector3.zero, 4 * rayDir, Color.red * (rayDir.magnitude < 1 ? 0.1f: 1.0f));
            
        }

        Debug.DrawRay(Vector3.zero, dir * 5, Color.green);
        Debug.DrawRay(Vector3.zero, newDir2 * 5, Color.magenta);
    }
}
