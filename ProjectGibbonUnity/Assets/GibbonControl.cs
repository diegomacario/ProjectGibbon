﻿using ImGuiNET;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Collections;

namespace Wolfire {
public class GibbonControl : MonoBehaviour
{
    public GameObject gibbon;
    public GameObject display_gibbon;
    public GameObject point_prefab;
    static public GameObject point_prefab_static;
    
    class HandState {
        public float3 pos;
        public bool gripping;
    }
    float last_x;

    HandState[] hands;
    int next_hand = 1;
    const float arm_length = 0.8f;
    
    float measured_arm_length;

    class VerletPoint {
        public float3 bind_pos;
        public float3 pos;
        public float3 old_pos;
        public float3 temp;
        public float mass;
        public bool pinned;
        public string name;
        public Transform widget;
    }

    class VerletBone {
        public int2 points;
        public float2 length;
        public string name;
        public bool enabled;
    }

    class BindPart {
        public Transform transform;
        public float4x4 bind_mat;
        public quaternion bind_rot;
        public float3 bind_pos;
    }

    class BindParts {
        public BindPart chest = new BindPart();
        public BindPart arm_top_l = new BindPart();
        public BindPart arm_bottom_l = new BindPart();
        public BindPart arm_top_r = new BindPart();
        public BindPart arm_bottom_r = new BindPart();
        public BindPart head = new BindPart();
        public BindPart belly = new BindPart();
        public BindPart pelvis = new BindPart();
        public BindPart leg_top_l = new BindPart();
        public BindPart leg_bottom_l = new BindPart();
        public BindPart leg_top_r = new BindPart();
        public BindPart leg_bottom_r = new BindPart();
    }

    BindParts bind_parts = new BindParts();

    class VerletSystem {
        public List<VerletPoint> points = new List<VerletPoint>();
        public List<VerletBone> bones = new List<VerletBone>();
        
        public void AddPoint(float3 pos, string name){
            var point = new VerletPoint();
            point.bind_pos = pos;
            point.pos = pos;
            point.old_pos = pos;
            point.pinned = false;
            point.name = name;
            point.mass = 1.0f;
            point.widget = Instantiate(point_prefab_static, point.pos, Quaternion.identity).transform;
            points.Add(point);
        }
    
        public void AddBone(string name, int a, int b) {
            var bone = new VerletBone();
            bone.points[0] = a;
            bone.points[1] = b;
            bone.length = math.distance(points[b].pos, points[a].pos);
            bone.name = name;
            bone.enabled = true;
            bones.Add(bone);
        }

        public void DrawBones(Color color){
            for(int i=0; i<bones.Count; ++i){
                if(bones[i].enabled){
                    DebugDraw.Line(points[bones[i].points[0]].pos,points[bones[i].points[1]].pos, color, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                }
            }       
        }

        public void StartSim(float step){
            float time_sqrd = step * step;
            var acc = new float3(0, Physics.gravity[1], 0);
            for(int i=0; i<points.Count; ++i){
                var point = points[i];
                point.temp = point.pos;
                if(!point.pinned){
                    point.pos = point.pos + (point.pos - point.old_pos) +  acc * time_sqrd;
                }
            }
        }

        public void Constraints() {
            for(int i=0; i<bones.Count; ++i){
                if(!bones[i].enabled){
                    continue;
                }
                var bone = bones[i];
                int num_pinned = 0;
                if(points[bone.points[0]].pinned) {++num_pinned;}
                if(points[bone.points[1]].pinned) {++num_pinned;}
                if(num_pinned < 2){
                    float curr_len = math.distance(points[bone.points[0]].pos, points[bone.points[1]].pos);
                    if(curr_len != 0f){
                        if(num_pinned == 1){
                            int pinned = bone.points[1];
                            int unpinned = bone.points[0];
                            if(points[bone.points[0]].pinned){
                                pinned = bone.points[0];
                                unpinned = bone.points[1];
                            }
                            var unpinned_point = points[unpinned];
                            if(curr_len < bone.length[0]){
                                unpinned_point.pos = points[pinned].pos + (points[unpinned].pos - points[pinned].pos) / curr_len * bone.length[0];
                            } else if(curr_len > bone.length[1]){
                                unpinned_point.pos = points[pinned].pos + (points[unpinned].pos - points[pinned].pos) / curr_len * bone.length[1];
                            }
                        } else {
                            var offset = (points[bone.points[1]].pos - points[bone.points[0]].pos) / curr_len;
                            float rel_mass = points[bone.points[1]].mass / (points[bone.points[0]].mass + points[bone.points[1]].mass);
                            var mid = points[bone.points[0]].pos * (1f-rel_mass) + points[bone.points[1]].pos * rel_mass;
                            if(curr_len < bone.length[0]){
                                points[bone.points[0]].pos = mid - offset * bone.length[0] * rel_mass;
                                points[bone.points[1]].pos = mid + offset * bone.length[0] * (1f-rel_mass);
                            } else if(curr_len > bone.length[1]){
                                points[bone.points[0]].pos = mid - offset * bone.length[1] * rel_mass;
                                points[bone.points[1]].pos = mid + offset * bone.length[1] * (1f-rel_mass);
                            }
                        }
                    }
                }
            }
        }

        public void EndSim() {
            for(int i=0; i < points.Count; ++i) {
                var point = points[i];
                point.old_pos = point.temp;
            }
        }

        public void Step(float step) {
            StartSim(step);
            Constraints();
            EndSim();
        }
    }

    VerletSystem pendulum = new VerletSystem();
    VerletSystem arms = new VerletSystem();
    VerletSystem complete = new VerletSystem();

    void SetBindPart(BindPart part, Transform transform){
        part.transform = transform;
        part.bind_mat = transform.localToWorldMatrix;
        part.bind_pos = transform.position;
        part.bind_rot = transform.rotation;
    }

    // Start is called before the first frame update
    
    class TwoBoneIK {
        public float3[] points;
        public TwoBoneIK() {
            points = new float3[3];
        }
    }

    TwoBoneIK arm_ik = new TwoBoneIK();
    TwoBoneIK leg_ik = new TwoBoneIK();

    float3 simple_pos;
    float3 target_com;
    float3 simple_vel = float3.zero;
    void Start() {
        point_prefab_static = point_prefab;

        simple_pos = gibbon.transform.position;
        simple_pos[1] = 0f;
        simple_pos[2] = 0f;
        hands = new HandState[2];
        for(int i=0; i<2; ++i){
            hands[i] = new HandState();
            hands[i].pos = simple_pos;
            hands[i].gripping = true;
        }
        hands[1].gripping = false;

        pendulum_center = simple_pos;
        //pendulum_length = 
        
        //gibbon.GetComponent<Animator>().runtimeAnimatorController.animationClips[3].SampleAnimation(gibbon, 0.0f);
        //character.GetTransforms(gibbon.transform.Find("rig/root"));
        //character.Draw();            
        
        var root = GameObject.Find("points").transform;
        var neck = root.Find("neck");
        var stomach = root.Find("stomach");
        var pelvis = root.Find("pelvis");
        var groin = root.Find("groin");
        var head = root.Find("head");
        var shoulder = root.Find("shoulder");
        var elbow = root.Find("elbow");
        var grip = root.Find("grip");
        var hip = root.Find("hip");
        var knee = root.Find("knee");
        var foot = root.Find("foot");        

        //DebugDraw.Line(neck.position, shoulder.position, Color.white, DebugDraw.Lifetime.Persistent, DebugDraw.Type.Xray);
        //DebugDraw.Line(shoulder.position, elbow.position, Color.white, DebugDraw.Lifetime.Persistent, DebugDraw.Type.Xray);
        //DebugDraw.Line(elbow.position, grip.position, Color.white, DebugDraw.Lifetime.Persistent, DebugDraw.Type.Xray);
        
        SetBindPart(bind_parts.head, display_gibbon.transform.Find("DEF-head"));
        SetBindPart(bind_parts.chest, display_gibbon.transform.Find("DEF-chest"));
        SetBindPart(bind_parts.belly, display_gibbon.transform.Find("DEF-belly"));
        SetBindPart(bind_parts.pelvis, display_gibbon.transform.Find("DEF-pelvis"));
        SetBindPart(bind_parts.arm_top_l, display_gibbon.transform.Find("DEF-upper_arm_L"));
        SetBindPart(bind_parts.arm_bottom_l, display_gibbon.transform.Find("DEF-forearm_L"));
        SetBindPart(bind_parts.arm_top_r, display_gibbon.transform.Find("DEF-upper_arm_R"));
        SetBindPart(bind_parts.arm_bottom_r, display_gibbon.transform.Find("DEF-forearm_R"));
        SetBindPart(bind_parts.leg_top_l, display_gibbon.transform.Find("DEF-thigh_L"));
        SetBindPart(bind_parts.leg_bottom_l, display_gibbon.transform.Find("DEF-shin_L"));
        SetBindPart(bind_parts.leg_top_r, display_gibbon.transform.Find("DEF-thigh_R"));
        SetBindPart(bind_parts.leg_bottom_r, display_gibbon.transform.Find("DEF-shin_R"));

        elbow.position = bind_parts.arm_bottom_r.transform.position;

        measured_arm_length = Vector3.Distance(shoulder.position, elbow.position) + Vector3.Distance(elbow.position, grip.position);
                  
        arm_ik.points[0] = shoulder.position;
        arm_ik.points[1] = elbow.position;
        arm_ik.points[2] = grip.position;
        
        leg_ik.points[0] = hip.position;
        leg_ik.points[1] = bind_parts.leg_bottom_r.transform.position;
        leg_ik.points[2] = foot.position;

        arms.AddPoint(shoulder.position, "shoulder_r");
        arms.AddPoint(grip.position, "hand_r");
        arms.AddPoint((shoulder.position+Vector3.right * (neck.position[0] - shoulder.position[0])*2f), "shoulder_l");
        arms.AddPoint((grip.position+Vector3.right * (neck.position[0] - grip.position[0])*2f), "hand_l");
        arms.AddPoint(new float3(neck.position[0], hip.position[1], neck.position[2]), "body");
        arms.points[0].mass = 2f;
        arms.points[2].mass = 2f;
        arms.points[4].mass = 4f;
        
        arms.AddBone("arm_r", 0, 1);
        arms.bones[arms.bones.Count-1].length[1] = measured_arm_length;
        arms.bones[arms.bones.Count-1].length[0] *= 0.4f; // Allow arm to flex
        arms.AddBone("arm_l", 2, 3);
        arms.bones[arms.bones.Count-1].length[1] = measured_arm_length;
        arms.bones[arms.bones.Count-1].length[0] *= 0.4f;
        arms.AddBone("tri_top", 0, 2);
        arms.AddBone("tri_r", 0, 4);
        arms.AddBone("tri_l", 2, 4);
        
        complete.AddPoint(shoulder.position, "shoulder_r");
        complete.AddPoint(grip.position, "hand_r");
        complete.AddPoint((shoulder.position+Vector3.right * (neck.position[0] - shoulder.position[0])*2f), "shoulder_l");
        complete.AddPoint((grip.position+Vector3.right * (neck.position[0] - grip.position[0])*2f), "hand_l");
        complete.AddPoint(new float3(neck.position[0], hip.position[1], neck.position[2]), "body");
        complete.AddPoint(head.position, "head");
        complete.AddPoint(neck.position, "neck");
        complete.AddPoint(stomach.position, "stomach");
        complete.AddPoint(pelvis.position, "hip");
        complete.AddPoint(groin.position, "groin");
        complete.AddPoint(hip.position, "hip_r");
        complete.AddPoint(foot.position, "foot_r");
        complete.AddPoint(hip.position+Vector3.right * (neck.position[0] - hip.position[0])*2f, "hip_l");
        complete.AddPoint(foot.position+Vector3.right * (neck.position[0] - foot.position[0])*2f, "foot_l");
        
        complete.AddBone("arm_r", 0, 1);
        complete.bones[complete.bones.Count-1].length[1] = measured_arm_length;
        complete.bones[complete.bones.Count-1].length[0] *= 0.4f; // Allow arm to flex
        complete.AddBone("arm_l", 2, 3);
        complete.bones[complete.bones.Count-1].length[1] = measured_arm_length;
        complete.bones[complete.bones.Count-1].length[0] *= 0.4f;
        complete.AddBone("head", 5, 6);
        complete.AddBone("chest", 6, 7);
        complete.AddBone("belly", 7, 8);
        complete.AddBone("pelvis", 8, 9);
        complete.AddBone("leg_r", 10, 11);
        complete.bones[complete.bones.Count-1].length[0] *= 0.4f;
        complete.AddBone("leg_l", 12, 13);
        complete.bones[complete.bones.Count-1].length[0] *= 0.4f;


        pendulum.AddPoint(new float3(pendulum_center), "pendulum_axis");
        pendulum.AddPoint(new float3(pendulum_center + new float3(0, -pendulum_length, 0)), "pendulum_end");
        pendulum.AddPoint(new float3(pendulum_center), "pendulum_next_grip");
        pendulum.points[0].pinned = true;
        pendulum.points[2].pinned = true;
        pendulum.AddBone("pendulum", 0, 1);
        pendulum.AddBone("pendulum_next", 2, 1);
        
    }
    

        static void DrawCircle(Vector3 pos, float radius){
        int num_segments = 32;
        for(int i=1; i<num_segments+1; ++i){
            float interp = (i-1)/(float)num_segments * Mathf.PI * 2f;
            float interp2 = i/(float)num_segments * Mathf.PI * 2f;
            var temp_pos = pos + (Vector3.right * Mathf.Sin(interp) - Vector3.up * Mathf.Cos(interp))*radius;
            var temp_pos2 = pos + (Vector3.right * Mathf.Sin(interp2) - Vector3.up * Mathf.Cos(interp2))*radius;
            DebugDraw.Line(temp_pos, temp_pos2, new Color(1.0f, 0.0f, 0.0f, 1.0f), DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Normal);
        }
    }

    Vector3 clicked_point;
    float3 com;

    float3 pendulum_center;
    float pendulum_length = 0.9f;

    void GripPoint(float3 pos){
        hands[next_hand].pos = pos;
        hands[next_hand].gripping = true;
        hands[1-next_hand].gripping = false;
        next_hand = 1 - next_hand;
        pendulum.points[0].pos = pos;
        pendulum.points[2].pos = pos;
        pendulum.bones[0].enabled = true;
        pendulum.bones[1].enabled = true;
        float len = math.distance(pendulum.points[0].pos, pendulum.points[1].pos);
        pendulum.bones[0].length[0] = len;
        pendulum.bones[1].length[0] = len;
        pendulum.bones[0].length[1] = len;
        pendulum.bones[1].length[1] = len;
    }
    
    float swing_time = 0f;

    static float GetAngleGivenSides(float a, float b, float c){
        // law of cosines:
        // c*c = a*a + b*b - 2*a*b*cos(C)
        // c*c - a*a - b*b = -2*a*b*cos(C)
        // (c*c - a*a - b*b) / (-2*a*b) = cos(C)
        // C = acos((c*c - a*a - b*b) / (-2*a*b))
        var top = (c*c - a*a - b*b);
        var divisor = (-2*a*b);
        if(divisor==0f){
            return 0f;
        }
        return math.acos(math.clamp(top / divisor, -1f, 1f));
    }

    void ApplyBound(BindPart part, float3 forward, float3 bind_forward, int start, int end){
        var up = math.normalize(complete.points[end].pos  - complete.points[start].pos);
        var bind_up = math.normalize(complete.points[end].bind_pos  - complete.points[start].bind_pos);       
        var mid = (complete.points[end].pos + complete.points[start].pos)/2.0f;
        var bind_mid = (complete.points[end].bind_pos + complete.points[start].bind_pos)/2.0f;
        
        var rotation = Quaternion.LookRotation(up, forward) * 
                       Quaternion.Inverse(Quaternion.LookRotation(bind_up, bind_forward));
        part.transform.rotation = rotation * part.bind_rot;
        part.transform.position = mid + (float3)(rotation * (part.bind_pos - bind_mid));
    }

    static void ApplyTwoBoneIK(int start, int end, float3 forward, TwoBoneIK ik, BindPart top, BindPart bottom, List<VerletPoint> points, float3 old_axis, float3 axis){
            var shoulder_offset = (points[start].pos - points[start].bind_pos);
            var shoulder_rotation = Quaternion.LookRotation(points[end].pos - points[start].pos, forward) * Quaternion.Inverse(Quaternion.LookRotation(points[end].bind_pos - points[start].bind_pos, Vector3.forward));
        
            float dist_a = math.distance(ik.points[0], ik.points[1]);
            float dist_b = math.distance(ik.points[1], ik.points[2]);
            float dist_c = math.distance(points[start].pos, points[end].pos);
            float old_dist_c = math.distance(ik.points[0], ik.points[2]);
            var old_elbow_angle = GetAngleGivenSides(dist_a, dist_b, old_dist_c);
            var old_shoulder_angle = GetAngleGivenSides(old_dist_c, dist_a, dist_b);
            var elbow_angle = GetAngleGivenSides(dist_a, dist_b, dist_c);
            var shoulder_angle = GetAngleGivenSides(dist_c, dist_a, dist_b);
            DebugText.AddVar("elbow_angle", elbow_angle, 0.5f);
            DebugText.AddVar("shoulder_angle", shoulder_angle, 0.5f);

            // Elbow axis is perpendicular to arm direction and vector from middle of arm to base of neck
            shoulder_rotation = Quaternion.AngleAxis(shoulder_angle * Mathf.Rad2Deg, axis) * shoulder_rotation * 
                                Quaternion.Inverse(Quaternion.AngleAxis(old_shoulder_angle * Mathf.Rad2Deg, old_axis));
            
            top.transform.position = top.bind_pos + shoulder_offset;
            top.transform.rotation = shoulder_rotation * top.bind_rot;
        
            var elbow = top.transform.position + top.transform.rotation * Quaternion.Inverse(top.bind_rot) * (bottom.bind_pos - top.bind_pos);
            bottom.transform.position = elbow;
            bottom.transform.rotation = Quaternion.AngleAxis(elbow_angle * Mathf.Rad2Deg, axis) * shoulder_rotation * 
                                        Quaternion.Inverse(Quaternion.AngleAxis(old_elbow_angle * Mathf.Rad2Deg, old_axis)) * bottom.bind_rot;
        
    }

    // Update is called once per frame
    void Update()
    {
        // Gibbon top speed brachiation is about 15 m/s
        // Can leap up to 8 meters
        // Gibbon wrist is ball and socket
        // Weigh ~7kg, 90 cm tall
        // Gibbons much slower on legs than trees
        // Leg jump speed 8.3 m/s (most force from swinging arms)
        // Ground speed up to 4 m/s?
        // Continous contact can go up to 4 m/s, ricochetal can go as low as 2.5 m/s
        // Ricochetal contact time ~0.5s at 3 m/s, 0.25s at 6 m/s
        // Continuous contact optimized when spaced slightly closer than full arm spread of animal, so around 1.2 m
        // Usually has margin of error, swings with arm not fully extended at high speed
        // Pendulum from 80 degrees w length 1 m has base speed of 4 m/s
                
        /*
        var editor = GetComponent<AnimationEditor>();
        foreach(var pose in editor.poses){
            if(pose.name == "Hang_pole"){
                AnimationEditor.ApplyPose(gibbon.transform, pose);                
            }
        }*/
        
        const bool manually_drag_points = false;
        if(manually_drag_points){
            for(int i=0; i<arms.points.Count; ++i){
                arms.points[i].pos = arms.points[i].widget.position;
            }
            arms.Constraints();
            for(int i=0; i<arms.points.Count; ++i){
                arms.points[i].widget.position = arms.points[i].pos;
            }
        }

        const bool map_complete_to_arms = true;
        if(map_complete_to_arms){
            var bind_mid = (arms.points[0].bind_pos + arms.points[2].bind_pos + arms.points[4].bind_pos)/3.0f;
            var mid = (arms.points[0].pos + arms.points[2].pos + arms.points[4].pos)/3.0f;
            var forward = math.normalize(math.cross(arms.points[0].pos - arms.points[2].pos, arms.points[0].pos - arms.points[4].pos));
            var bind_forward = math.normalize(math.cross(arms.points[0].bind_pos - arms.points[2].bind_pos, arms.points[0].bind_pos - arms.points[4].bind_pos));
            var up = math.normalize((arms.points[0].pos + arms.points[2].pos)/2.0f - arms.points[4].pos);
            var bind_up = math.normalize((arms.points[0].bind_pos + arms.points[2].bind_pos)/2.0f - arms.points[4].bind_pos);
        
            complete.points[0].pos = arms.points[0].pos;
            complete.points[1].pos = arms.points[1].pos;
            complete.points[2].pos = arms.points[2].pos;
            complete.points[3].pos = arms.points[3].pos;
            complete.points[0].pinned = true;
            complete.points[1].pinned = true;
            complete.points[2].pinned = true;
            complete.points[3].pinned = true;

            var chest_rotation = math.mul(quaternion.LookRotation(forward, up), 
                                          math.inverse(quaternion.LookRotation(bind_forward, bind_up)));

            for(int i=5; i<14; ++i){
                complete.points[i].pos = mid + math.mul(chest_rotation, (complete.points[i].bind_pos - bind_mid));
                complete.points[i].pinned = true;
            }

            // Leg motion
            for(int i=0; i<2; ++i){
                complete.points[11+i*2].pos += up * (0.35f + 0.15f * math.sin((swing_time + i*1.0f)*math.PI));
            }

            complete.Constraints();
        }

        {
            var bind_mid = (complete.points[0].bind_pos + complete.points[2].bind_pos + complete.points[9].bind_pos)/3.0f;
            var mid = (complete.points[0].pos + complete.points[2].pos + complete.points[9].pos)/3.0f;
            var forward = -math.normalize(math.cross(complete.points[0].pos - complete.points[2].pos, complete.points[0].pos - complete.points[9].pos));
            var bind_forward = -math.normalize(math.cross(complete.points[0].bind_pos - complete.points[2].bind_pos, complete.points[0].bind_pos - complete.points[9].bind_pos));
            var up = math.normalize((complete.points[0].pos + complete.points[2].pos)/2.0f - complete.points[9].pos);
            var bind_up = math.normalize((complete.points[0].bind_pos + complete.points[2].bind_pos)/2.0f - complete.points[9].bind_pos);
        
            //DebugDraw.Line(mid, mid+forward, Color.blue, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
            //DebugDraw.Line(mid, mid+up, Color.green, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);

            bind_parts.chest.transform.rotation = Quaternion.LookRotation(forward, up) * 
                                                  Quaternion.Inverse(Quaternion.LookRotation(bind_forward, bind_up)) * 
                                                  bind_parts.chest.bind_rot;
            bind_parts.chest.transform.position = mid + 
                                                  (float3)(Quaternion.LookRotation(forward, up) * 
                                                  Quaternion.Inverse(Quaternion.LookRotation(bind_forward, bind_up)) * 
                                                  (bind_parts.chest.bind_pos - bind_mid));
                                            
            ApplyBound(bind_parts.head, forward, bind_forward, 5, 6);
            ApplyBound(bind_parts.chest, forward, bind_forward, 6, 7);
            ApplyBound(bind_parts.belly, forward, bind_forward, 7, 8);
            ApplyBound(bind_parts.pelvis, forward, bind_forward, 8, 9);
            ApplyBound(bind_parts.leg_top_r, forward, bind_forward, 10, 11);
            ApplyBound(bind_parts.leg_top_l, forward, bind_forward, 12, 13);

            for(int i=0; i<2; ++i){
                var top = bind_parts.arm_top_r;
                var bottom = bind_parts.arm_bottom_r;
                if(i==1){
                    top = bind_parts.arm_top_l;
                    bottom = bind_parts.arm_bottom_l;
                }

                var points = complete.points;
                int start = i*2;
                int end = i*2+1;
                var old_axis = math.normalize(math.cross((points[end].bind_pos+points[start].bind_pos)*0.5f - (points[2].bind_pos + points[0].bind_pos) * 0.5f, points[start].bind_pos - points[end].bind_pos));//shoulder_rotation * Vector3.forward;// math.normalize(shoulder_rotation * math.cross(left_arm.points[2] - left_arm.points[1], left_arm.points[1] - left_arm.points[0]));
                var axis = math.normalize(math.cross((points[end].pos+points[start].pos)*0.5f - (points[2].pos + points[0].pos) * 0.5f, points[start].pos - points[end].pos));//shoulder_rotation * Vector3.forward;// math.normalize(shoulder_rotation * math.cross(left_arm.points[2] - left_arm.points[1], left_arm.points[1] - left_arm.points[0]));
            
                ApplyTwoBoneIK(start, end, forward, arm_ik, top, bottom, complete.points, old_axis, axis);
            }

            for(int i=0; i<2; ++i){
                var top = bind_parts.leg_top_r;
                var bottom = bind_parts.leg_bottom_r;
                if(i==1){
                    top = bind_parts.leg_top_l;
                    bottom = bind_parts.leg_bottom_l;
                }
            
                var points = complete.points;
                int start = i*2+10;
                int end = i*2+1+10;

                var bind_shoulder_rotation = Quaternion.LookRotation(points[end].bind_pos - points[start].bind_pos, Vector3.forward);
                var shoulder_rotation = Quaternion.LookRotation(points[end].pos - points[start].pos, forward) * bind_shoulder_rotation;
        
                var old_axis = bind_shoulder_rotation * Vector3.right;
                var axis = shoulder_rotation * Vector3.right;

                ApplyTwoBoneIK(start, end, forward, leg_ik, top, bottom, complete.points, old_axis, axis);
            }
        }

        //gibbon.transform.position = mid;
        //gibbon.transform.rotation = Quaternion.LookRotation(forward, up);
        
        //arms.DrawBones(Color.white);
        //pendulum.DrawBones(Color.blue);
        //complete.DrawBones(Color.green);
        
        if(!hands[0].gripping && !hands[1].gripping){
            var vel = (pendulum.points[1].pos - pendulum.points[1].old_pos) / (Time.fixedDeltaTime * 0.1f);
            var perp = new float3(vel[1], -vel[0], 0f);
            if(perp[1] < 0){
                perp = -perp;
            }
            perp = math.normalize(perp);
            float len = -pendulum.points[1].pos[1] / perp[1];
            var point = pendulum.points[1].pos + perp * len;
            DebugDraw.Line(pendulum.points[1].pos, point, Color.red, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
        }

        if(Input.GetKeyDown(KeyCode.Space)){
            if(hands[0].gripping || hands[1].gripping ){
                hands[0].gripping = false;
                hands[1].gripping = false;
                pendulum.bones[0].enabled = false;
                pendulum.bones[1].enabled = false;
            } else {
                var vel = (pendulum.points[1].pos - pendulum.points[1].old_pos) / (Time.fixedDeltaTime * 0.1f);
                var perp = new float3(vel[1], -vel[0], 0f);
                if(perp[1] < 0){
                    perp = -perp;
                }
                perp = math.normalize(perp);
                float len = -pendulum.points[1].pos[1] / perp[1];
                /*if(len > arm_length){
                    len = arm_length;
                }*/
                var point = pendulum.points[1].pos + perp * len;
                point[1] = 0f;
                GripPoint(point);
            }
        }

        float total_mass = 0f;
        var old_com = com;
        com = float3.zero;
        for(int i=0; i<arms.points.Count; ++i){
            com += arms.points[i].pos * arms.points[i].mass;
            total_mass += arms.points[i].mass;
        }
        com /= total_mass;
        //DebugDraw.Sphere(com, Color.green, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray );
                
        if(ImGui.Begin("Gibbon")){
            //ImGui.Text($"pendulum_length: {math.distance(arms.points[1].pos, com)}");
            //DebugDraw.Line(arms.points[1].pos, com, Color.green, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray );    
            float3 vel = (pendulum.points[1].pos - pendulum.points[1].old_pos) / (Time.fixedDeltaTime * 0.1f);
            float kinetic_energy = math.lengthsq(vel) / 2.0f;
            ImGui.Text($"kinetic energy: {kinetic_energy}");   
            float height = pendulum.points[1].pos[1] - (pendulum.points[0].pos[1] - pendulum.bones[0].length[1]);
            float potential_energy = height*-Physics.gravity[1];
            ImGui.Text($"potential energy: {potential_energy}");   
            ImGui.Text($"total energy: {kinetic_energy + potential_energy}");  
            ImGui.Text($"horz speed: {math.abs(vel[0])}");  
            ImGui.Text($"horz speed: {math.abs(simple_vel[0])}");  
            ImGui.Text($"swing time: {swing_time}");  
        }
        ImGui.End();

        if(Input.GetKeyDown(KeyCode.Tab)){
            Time.timeScale = (Time.timeScale == 1.0f)?0.1f:1.0f;
        }

        const bool pendulum_test = false;
        if(pendulum_test){
            float pendulum_period = math.PI * 2.0f * math.sqrt(pendulum_length / -Physics.gravity[1]);
            float pendulum_angle = math.sin(math.PI * 2.0f / pendulum_period * Time.time);
            var pendulum_pos = new float3(math.sin(pendulum_angle), -math.cos(pendulum_angle), 0f) * pendulum_length + pendulum_center;
            DebugDraw.Sphere(pendulum_pos, Color.blue, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray );
        
            DebugDraw.Line(pendulum_center, pendulum_pos, Color.blue, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray );    
        }

        //DebugDraw.Sphere(pendulum.points[1].pos, Color.blue, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray );
        
    }
    
    float grabbed_time = 0f;

    float3 MoveTowards(float3 a, float3 b, float max_dist){
        float len = math.distance(a,b);
        if(len < max_dist){
            return b;
        } else {
            return a + (b-a)/len*max_dist;
        }
    }

    void Step(float step) {
        float horz_input = 0f;
        float vert_input = 0f;
        if(Input.GetKey(KeyCode.D)){
            horz_input = 1f;
        }
        if(Input.GetKey(KeyCode.A)){
            horz_input = -1f;
        }
        if(Input.GetKey(KeyCode.W)){
            vert_input = 1f;
        }
        if(Input.GetKey(KeyCode.S)){
            vert_input = -1f;
        }

        
        var old_pos = simple_pos;
        simple_vel[0] += horz_input * Time.deltaTime * 5f;
        simple_vel[0] = math.clamp(simple_vel[0], -10f, 10f);
        simple_pos += simple_vel * Time.deltaTime;
        // Adjust amplitude and time scale based on speed
        float amplitude = math.pow(math.abs(simple_vel[0])/10f + 1f, 0.8f)-1f+0.1f;
        float min_height = -1f + amplitude * 0.25f + math.max(0.0f, 0.1f - math.abs(simple_vel[0]) * 0.1f);
        var old_swing_time = swing_time;
        float swing_speed_mult = 8f/(math.PI*2f);
        swing_time += step*swing_speed_mult;
        if(math.ceil(old_swing_time) != math.ceil(swing_time)){
            next_hand = 1-next_hand;
        }
        simple_pos[1] = (min_height + (math.sin((swing_time-0.1f) * (math.PI*2f))+1f)*amplitude) * pendulum_length;
        target_com = simple_pos - simple_vel * 0.05f;
        target_com[0] += (math.cos((swing_time-0.1f) * (math.PI*2f)))* pendulum_length * 0.5f * math.clamp(simple_vel[0] * 0.5f, -1f, 1f) * math.max(0f, 1f - math.abs(simple_vel[0])*2f);
        //DebugDraw.Sphere(simple_pos, Color.green, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray );
        //DebugDraw.Line(old_pos, simple_pos, Color.green, DebugDraw.Lifetime.Persistent, DebugDraw.Type.Xray );    
        float next_trough_time = ((math.ceil(swing_time)-0.25f))/swing_speed_mult;
        hands[next_hand].pos = simple_pos + simple_vel * (next_trough_time-Time.time);
        hands[next_hand].pos[1] = 0.0f;
        DebugDraw.Sphere(hands[0].pos, Color.green, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray );
        DebugDraw.Sphere(hands[1].pos, Color.blue, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray );
        
        
        // If movement key is held, then look at trajectory
        // If we have a good jump trajectory, then let go and follow trajectory and grab at best grab point
        // Otherwise if we're at the end of our swing and can reach the next handhold, then start a new swing (going backwards slightly for momentum)
        // Otherwise swing more to get momentum
        
        bool check_for_jump = false;
        if(check_for_jump){
            if(horz_input != 0.0f && (hands[0].gripping || hands[1].gripping)) {
                var temp_pos = pendulum.points[1].pos;
                float3 vel = (pendulum.points[1].pos - pendulum.points[1].old_pos) / (Time.fixedDeltaTime * 0.1f);
                for (int i = 0; i < 30; ++i) {
                    temp_pos += vel * Time.fixedDeltaTime;
                    DebugDraw.Sphere(temp_pos, Color.red, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);

                    var perp = new float3(vel[1], -vel[0], 0f);
                    if (perp[1] < 0) {
                        perp = -perp;
                    }
                    bool any_good_trajectory = false;
                    if (perp[1] != 0) {
                        perp = math.normalize(perp);
                        float len = -temp_pos[1] / perp[1];
                        var point = temp_pos + perp * len;
                        var color = Color.red;
                        float dist = math.distance(point, temp_pos);
                        // Highlight grip if it will be within reach and will be going in the same direction
                        if (dist < 1.0f && dist > 0.75f && perp[0] * vel[0] > 0.0f && math.abs(point[0] - hands[1-next_hand].pos[0]) > arm_length * 2f) {
                            color = Color.green;
                            any_good_trajectory = true;
                        }
                        DebugDraw.Line(temp_pos, point, color, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
                        if(any_good_trajectory){
                            // let go
                            hands[0].gripping = false;
                            hands[1].gripping = false;
                            pendulum.bones[0].enabled = false;
                            pendulum.bones[1].enabled = false;
                            Debug.Log("LET GO");
                            break;
                        }
                    }

                    vel += (float3)Physics.gravity * 0.1f;
                }
            }
        }

        bool use_pendulum = false;
        if(use_pendulum){
            for (int i=0; i<10; ++i){
                if(!hands[0].gripping && !hands[1].gripping){
                    float3 vel = (pendulum.points[1].pos - pendulum.points[1].old_pos) / (Time.fixedDeltaTime * 0.1f);
                    var temp_pos = pendulum.points[1].pos;

                    var perp = new float3(vel[1], -vel[0], 0f);
                    if (perp[1] < 0) {
                        perp = -perp;
                    }
                    if (perp[1] != 0) {
                        perp = math.normalize(perp);
                        float len = -temp_pos[1] / perp[1];
                        var point = temp_pos + perp * len;
                        var color = Color.red;
                        float dist = math.distance(point, temp_pos);
                        // Highlight grip if it will be within reach and will be going in the same direction
                        if (dist < 1.0f && dist > 0.75f && perp[0] * vel[0] > 0.0f && math.abs(point[0] - hands[1-next_hand].pos[0]) > arm_length * 2f) {
                            color = Color.green;
                            GripPoint(point);
                            Debug.Log("GRAB");
                            break;
                        }
                        DebugDraw.Line(temp_pos, point, color, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
                    }
                }

            float temp_step = step * 0.1f;
            float time_sqrd = temp_step * temp_step;
            pendulum.StartSim(temp_step);
            pendulum.points[1].pos[0] += horz_input * time_sqrd * 2f;
            pendulum.Constraints();
            pendulum.EndSim();

            /*
            if(horz_input != 0.0f){
                var grab_point = (float3)hands[1-next_hand].pos + new float3(arm_length * 1.8f * horz_input,0,0);
                if(math.distance(pendulum.points[1].pos, grab_point) < arm_length){        
                    GripPoint(grab_point);
                }
            }*/
            
                for(int a=0; a<2; ++a){
                    for(int b=0; b<2; ++b){
                        pendulum.bones[a].length[b] = Mathf.MoveTowards(pendulum.bones[a].length[b], pendulum_length, temp_step);
                    }
                }
                pendulum.Constraints();

            /*
            if(vert_input != 0.0f){
                pendulum.bones[0].length[0] -= vert_input * temp_step;
                pendulum.bones[0].length[1] -= vert_input * temp_step;
                pendulum.bones[1].length[0] -= vert_input * temp_step;
                pendulum.bones[1].length[1] -= vert_input * temp_step;
                pendulum.Constraints();
            }*/

            const bool draw_com_path = false;
            if(draw_com_path){
                DebugDraw.Line(pendulum.points[1].old_pos, pendulum.points[1].pos, Color.green, DebugDraw.Lifetime.Persistent, DebugDraw.Type.Xray );
            }
        }
        }

        bool arms_map = true;
        if(arms_map){
            /*arms.points[1].pinned = hands[0].gripping;
            if(hands[0].gripping){
                arms.points[1].pos = pendulum.points[0].pos;
            } else {
                arms.points[1].pos = arms.points[0].pos + new float3(measured_arm_length,0,0);
            }
            arms.points[3].pinned = hands[1].gripping;
            if(hands[1].gripping){
                arms.points[3].pos = pendulum.points[2].pos;
            } else {
                arms.points[3].pos = arms.points[2].pos + new float3(measured_arm_length,0,0);
            }*/
            
            for(int i=0; i<2; ++i){
               // if((hands[i].pos[0] - simple_pos[0]) * simple_vel[0] > 0f){
                    arms.points[i*2+1].pos = MoveTowards(arms.points[i*2+1].pos, hands[i].pos, math.max(0f, math.cos((swing_time+0.35f+(1-i))*math.PI*1f)*0.5f+0.5f) * step * 5f);
                //}   
                arms.points[i*2+1].old_pos = math.lerp(arms.points[i*2+1].old_pos, arms.points[i*2+1].pos - simple_vel*step, 0.25f);
            }
            arms.StartSim(step);
            for(int j=0; j<4; ++j){
                float total_mass = 0f;
                var com = float3.zero;
                for(int i=0; i<arms.points.Count; ++i){
                    if(arms.points[i].pinned == false){
                        com += arms.points[i].pos * arms.points[i].mass;
                        total_mass += arms.points[i].mass;
                    }
                }
                com /= total_mass;
                var offset = (float3)target_com - com;
                /*if(arms.points[1].pinned && !arms.points[3].pinned){
                    arms.points[0].pos += offset * total_mass / arms.points[0].mass;
                } else if(arms.points[3].pinned && !arms.points[1].pinned){
                    arms.points[2].pos += offset * total_mass / arms.points[2].mass;
                } else {
                    var temp_offset = offset * total_mass / (arms.points[0].mass + arms.points[2].mass);
                    arms.points[0].pos += temp_offset;
                    arms.points[2].pos += temp_offset;
                }*/
                for(int i=0; i<arms.points.Count; ++i){
                    if(i!=1 && i!=3){
                        arms.points[i].pos += offset;
                    }
                }
                // Apply torque to keep torso upright
                float step_sqrd = step*step;
                float force = 20f;
                arms.points[4].pos[1] -= step_sqrd * force;
                arms.points[0].pos[1] += step_sqrd * force * 0.5f;
                arms.points[2].pos[1] += step_sqrd * force * 0.5f;
                arms.points[0].pos[2] -= step_sqrd * simple_vel[0] * 2.0f;
                arms.points[2].pos[2] += step_sqrd * simple_vel[0] * 2.0f;
                arms.points[4].pos[0] -= simple_vel[0] * step_sqrd * 2f; // Apply backwards force to maintain forwards tilt
                arms.Constraints();
            }
            arms.EndSim();
        }

        bool arms_sim = false;
        if(arms_sim){        
            float time_sqrd = step * step;
            if(hands[0].gripping){    
                //points[1].pos = hands[0].pos;
                arms.points[1].pinned = true;
            } else {
                arms.points[1].pinned = false;
            }
            if(hands[1].gripping){    
                //points[3].pos = hands[1].pos;
                arms.points[3].pinned = true;
            } else {
                arms.points[3].pinned = false;
            }
        
            arms.StartSim(step);

            const int swing_in_place = 0, continuous = 1, ricochetal = 2; 

            int state = continuous;
            if(horz_input != 0f){
                var swing_force = 3f;
                float offset = horz_input * time_sqrd * swing_force;
                switch(state){
                    case swing_in_place:
                        arms.points[0].pos[0] += offset;
                        arms.points[2].pos[0] += offset;
                        arms.points[4].pos[0] += offset;
                        break;
                    case continuous:
                        if(grabbed_time <= Time.time - 0.2f){
                            hands[next_hand].gripping = false;
                        }
                        arms.points[0].pos[0] += offset;
                        arms.points[2].pos[0] += offset;
                        arms.points[4].pos[0] += offset;
                        var grip_pos = (float3)hands[1-next_hand].pos;
                        grip_pos[0] += arm_length * horz_input * 1.8f;
                        
                        if(grabbed_time <= Time.time - 0.4f){
                            float max_grab_force = 20f; 
                            var offset_len = time_sqrd * max_grab_force;// math.min(1f, d_key_held*1f) * max_grab_force;
                            int hand_point = next_hand*2+1;
                            if(math.distance(arms.points[hand_point].pos, grip_pos) < 0.1f){
                                //points[hand_point].pos = grip_pos;
                                arms.points[hand_point].pinned = true;
                                hands[next_hand].gripping = true;
                                hands[next_hand].pos = grip_pos;
                                next_hand = 1-next_hand;
                                //hands[next_hand].gripping = false;
                                grabbed_time = Time.time;
                            } else {
                                arms.points[hand_point].pos += math.normalize((float3)grip_pos - arms.points[hand_point].pos) * offset_len;
                            }
                        }
                        break;
                    case ricochetal:
                        break;
                }
            } 
            arms.Constraints();
            arms.Constraints();
            arms.Constraints();
            arms.Constraints();
            arms.EndSim();
        }
        var cam_pos = Camera.main.transform.position;
        cam_pos[0] = simple_pos[0];
        Camera.main.transform.position = cam_pos;
    }

    private void FixedUpdate() {
        Step(Time.fixedDeltaTime);
    }
}
}