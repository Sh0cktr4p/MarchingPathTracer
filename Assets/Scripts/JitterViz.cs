using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class JitterViz : MonoBehaviour
{
    public Transform incidentTransform;

    // Start is called before the first frame update
    void Start()
    {
    }


    Vector3 JitterHemiSphere(Vector3 d, Vector3 i)
    {
        float theta = Random.value * Mathf.PI * 2;
        float phi = Random.value * Mathf.PI * 0.5f;

        float x = Mathf.Sin(theta) * Mathf.Sin(phi);
        float z = Mathf.Cos(theta) * Mathf.Sin(phi);
        float y = Mathf.Cos(phi);

        Vector3 rotAxis = new Vector3(d.z, 0, -d.x).normalized;

        Debug.DrawRay(transform.position, rotAxis * 1.5f, Color.red);

        float angle = Mathf.Acos(d.y);

        Quaternion q = Quaternion.AngleAxis(angle * 180 / Mathf.PI, rotAxis);

        Debug.Log("Axis: " + rotAxis + ", angle: " + angle);

        return q * new Vector3(x, y, z);
    }

    // Update is called once per frame
    void Update()
    {
        Vector3 incidentDir = incidentTransform.forward;

        Debug.DrawRay(transform.position, transform.up * 1.5f, Color.cyan);
        Debug.DrawRay(transform.position, incidentDir * 1.5f, Color.yellow);

        for (int i = 0; i < 100; i++)
        {
            Debug.DrawRay(transform.position, JitterHemiSphere(transform.up, incidentDir), Color.magenta, 0.5f);
        }
    }
}
