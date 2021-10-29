using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;
//using System.Diagnostics;
using UnityEngine;

public class PlayerMovement : MonoBehaviour
{
    private static string axis_leftHorizontal = "LHor";
    private static string axis_leftVertical = "LVer";
    private static string axis_Trigger = "Trig";
    private static string axis_rightHorizontal = "RHor";
    private static string axis_rightVertical = "RVer";
    private static string axis_crossHorizontal = "CHor";
    private static string axis_crossVertical = "CVer";
    
    private static string[] axes = {axis_leftHorizontal, axis_leftVertical, axis_Trigger, axis_rightHorizontal, axis_rightVertical, axis_crossHorizontal, axis_crossVertical};

    public float moveSpeed;
    public float turnSpeed;
    public float intensityGain;
    public float maxIntensity;
    public float apertureGain;
    public float focalDistGain;

    private bool modifyFeedback;

    private ImageEffect imageEffect;

    // Start is called before the first frame update
    void Start()
    {
        imageEffect = GetComponent<ImageEffect>();
        modifyFeedback = imageEffect.feedback;
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown("joystick button 4")){
            moveSpeed /= 2;
        }

        if (Input.GetKeyDown("joystick button 5"))
        {
            moveSpeed *= 2;
        }

        Vector3 translation = new Vector3(Input.GetAxis(axis_leftHorizontal), 0, Input.GetAxis(axis_leftVertical));
        Vector3 rotation = new Vector3(-Input.GetAxis(axis_rightVertical), Input.GetAxis(axis_rightHorizontal), Input.GetAxis(axis_Trigger));

        int apertureDelta = (int) Input.GetAxisRaw(axis_crossHorizontal);
        int focalDistDelta = (int) Input.GetAxisRaw(axis_crossVertical);

        imageEffect.aperture += apertureDelta * apertureGain;
        imageEffect.focalDist += focalDistDelta * focalDistGain;

        if(translation == Vector3.zero && rotation == Vector3.zero && apertureDelta == 0 && focalDistDelta == 0 && modifyFeedback){
            imageEffect.feedback = true;
        }
        else{
            imageEffect.feedback = false;
        }
        
        transform.Rotate(rotation * turnSpeed * Time.deltaTime, Space.Self);
        transform.Translate(translation * moveSpeed * Time.deltaTime, Space.Self);
        
        if(Input.GetKeyDown("joystick button 0"))
        {
            StartCoroutine(turnOnTheLights());
        }
        if(Input.GetKeyDown("joystick button 1")){
            //StopCoroutine(turnOnTheLights());
            imageEffect.intensity = 0;
        }
        if(Input.GetKeyDown("joystick button 2")){
            imageEffect.intensity *= 2;
        }
        if(Input.GetKeyDown("joystick button 3")){
            imageEffect.intensity /= 2;
        }
    }

    void TestInput(){
        foreach(string axis in axes){
            float val = Input.GetAxis(axis);

            if(val != 0){
                UnityEngine.Debug.Log(axis + ": " + val);
            }
        }
    }

    IEnumerator turnOnTheLights()
    {   
        while (imageEffect.intensity / maxIntensity < 0.99)
        {
            imageEffect.intensity = Mathf.Lerp(imageEffect.intensity, maxIntensity, intensityGain);
            yield return new WaitForEndOfFrame();
        }
    }
}
