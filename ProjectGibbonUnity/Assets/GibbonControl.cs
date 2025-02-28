﻿using ImGuiNET;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Collections;

namespace Wolfire {
public class GibbonControl : MonoBehaviour {
    public GameObject display_gibbon; // Character mesh and bone transforms
        
    // Skinning data
    class DisplayBone {
        public Transform transform;
        public quaternion bind_rot;
        public float3 bind_pos;

        public void Bind(Transform transform){
            this.transform = transform;
            bind_pos = transform.position;
            bind_rot = transform.rotation;
        }
    }

    class DisplayBody {
        public DisplayBone chest = new DisplayBone();
        public DisplayBone arm_top_l = new DisplayBone();
        public DisplayBone arm_bottom_l = new DisplayBone();
        public DisplayBone arm_top_r = new DisplayBone();
        public DisplayBone arm_bottom_r = new DisplayBone();
        public DisplayBone head = new DisplayBone();
        public DisplayBone belly = new DisplayBone();
        public DisplayBone pelvis = new DisplayBone();
        public DisplayBone leg_top_l = new DisplayBone();
        public DisplayBone leg_bottom_l = new DisplayBone();
        public DisplayBone leg_top_r = new DisplayBone();
        public DisplayBone leg_bottom_r = new DisplayBone();
    }

    DisplayBody display_body = new DisplayBody();
    
    // Particle simulation systems
    Verlet.System complete = new Verlet.System();
    Verlet.System branches = new Verlet.System();

    // Initial point positions for IK
    float3[] arm_ik = new float3[3];
    float3[] leg_ik = new float3[3];
    
    class MovementSystem {
        public float3 target_com;
        // 0 is right hand
        // 1 is left hand
        // 2 is right foot
        // 3 is left foot
        public float3[] limb_targets = new float3[4];
        public Verlet.System simple_rig = new Verlet.System();
        public float body_compress_amount = 0.0f;
    }

    MovementSystem walk = new MovementSystem();
    MovementSystem display = new MovementSystem();
    
    // Simple rig point ids
    // I just remembered the int values and used them directly, but this is good for reference
    const int p_shoulder_r = 0;
    const int p_hand_r =     1;
    const int p_shoulder_l = 2;
    const int p_hand_l =     3;
    const int p_base =       4;

    // Debug draw values and window
    class DebugInfo {
        public bool draw_walk_rig = false;
        public bool draw_ik_final = false;
        public bool draw_display_simple_rig = false;
        public bool draw_display_complete_rig = false;
        public bool draw_gibbon = true;
        public bool draw_elbow_ik_target = false;
        public bool draw_com_line = false;
        public bool draw_simple_point = false;
        int draw_layer = 0;
        string[] layers = new string[] {
            "Skinned",
            "IK",
            "Rig",
            "Simple Rig",
            "Source Rigs",
            "Particle"
        };
        public bool draw_hand_pull = false;
        public bool draw_trajectory = false;
        public bool draw_head_look = false;
        public bool draw_smoothing = false;
        public List<DebugDraw.DebugDrawLine> com_lines = new List<DebugDraw.DebugDrawLine>();

        public void DrawWindow() {
            if(ImGui.Begin("Debug Visualization")){
                if(ImGui.Combo("Draw", ref draw_layer, layers)){
                    draw_gibbon = (draw_layer==0);
                    draw_ik_final = (draw_layer==1);
                    draw_display_complete_rig = (draw_layer==2);
                    draw_display_simple_rig = (draw_layer==3);
                    draw_walk_rig = (draw_layer==4);
                    draw_simple_point = (draw_layer==5);
                }
                if(ImGui.Checkbox("Draw COM line", ref draw_com_line)){
                    if(!draw_com_line){
                        foreach(var line in com_lines){
                            DebugDraw.Remove(line);
                        }
                        com_lines.Clear();
                    }
                }
                ImGui.Checkbox("Draw path smoothing", ref draw_smoothing);
                ImGui.Checkbox("Draw hand pull", ref draw_hand_pull);
                ImGui.Checkbox("Draw head look", ref draw_head_look);
                ImGui.Checkbox("Draw elbow IK target", ref draw_elbow_ik_target);
                
                bool slow_motion = (Time.timeScale != 1.0f);
                if(ImGui.Checkbox("Slow motion [tab]", ref slow_motion)) {
                    Time.timeScale = (Time.timeScale == 1.0f)?0.1f:1.0f;
                }
            }
            ImGui.End();
        }
    }
    
    DebugInfo debug_info = new DebugInfo();
    
    // Simple character particle information
    float3 simple_pos;
    float3 simple_vel = float3.zero;

    // Time marker in animation sequences (so cycles can speed up and slow down without losing continuity)
    float walk_time = 0f;

    float3 look_target; // For head look IK

    float body_compress_amount = 0.0f; // Used to shorten the distance between neck and hips if necessary, e.g. during quadruped gallop

    // Various parameters that were being used to tune an animation
    float base_walk_height = 0.7f;
    float tilt_offset = 0.81f;
    float gallop_offset = 0.55f; // For biped gallop
    float quad_gallop_offset = 0.25f; // For quadruped gallop
    float quad_amount = 0.0f;

    void Start() {
        // Starting point
        simple_pos = display_gibbon.transform.position;
        simple_pos[1] = 0f;
        simple_pos[2] = 0f;

        // Init hand positions
        for(int i=0; i<4; ++i){
            display.limb_targets[i] = simple_pos;
            walk.limb_targets[i] = simple_pos;
        }

        // Get transforms of each skeleton point
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
        
        // Set up bind poses for each bone
        display_body.head.Bind(display_gibbon.transform.Find("DEF-head"));
        display_body.chest.Bind(display_gibbon.transform.Find("DEF-chest"));
        display_body.belly.Bind(display_gibbon.transform.Find("DEF-belly"));
        display_body.pelvis.Bind(display_gibbon.transform.Find("DEF-pelvis"));
        display_body.arm_top_l.Bind(display_gibbon.transform.Find("DEF-upper_arm_L"));
        display_body.arm_bottom_l.Bind(display_gibbon.transform.Find("DEF-forearm_L"));
        display_body.arm_top_r.Bind(display_gibbon.transform.Find("DEF-upper_arm_R"));
        display_body.arm_bottom_r.Bind(display_gibbon.transform.Find("DEF-forearm_R"));
        display_body.leg_top_l.Bind(display_gibbon.transform.Find("DEF-thigh_L"));
        display_body.leg_bottom_l.Bind(display_gibbon.transform.Find("DEF-shin_L"));
        display_body.leg_top_r.Bind(display_gibbon.transform.Find("DEF-thigh_R"));
        display_body.leg_bottom_r.Bind(display_gibbon.transform.Find("DEF-shin_R"));

        // Adjust elbow to match arm transform
        elbow.position = display_body.arm_bottom_r.transform.position;

        // Set up initial IK poses (just used to get bone lengths later)
        arm_ik[0] = shoulder.position;
        arm_ik[1] = elbow.position;
        arm_ik[2] = grip.position;
        
        leg_ik[0] = hip.position;
        leg_ik[1] = display_body.leg_bottom_r.transform.position;
        leg_ik[2] = foot.position;

        float measured_arm_length = Vector3.Distance(shoulder.position, elbow.position) + Vector3.Distance(elbow.position, grip.position);
            
        // Set up movement system particles and bones
        for(int i=0; i<2; ++i){
            Verlet.System new_simple_rig;
            switch(i){
                case 0:  new_simple_rig  = display.simple_rig; break;
                default:  new_simple_rig  = walk.simple_rig; break;
            }

            new_simple_rig.AddPoint(shoulder.position, "shoulder_r");
            new_simple_rig.AddPoint(grip.position, "hand_r");
            new_simple_rig.AddPoint((shoulder.position+Vector3.right * (neck.position[0] - shoulder.position[0])*2f), "shoulder_l");
            new_simple_rig.AddPoint((grip.position+Vector3.right * (neck.position[0] - grip.position[0])*2f), "hand_l");
            new_simple_rig.AddPoint(new float3(neck.position[0], hip.position[1], neck.position[2]), "body");
            new_simple_rig.points[0].mass = 2f;
            new_simple_rig.points[2].mass = 2f;
            new_simple_rig.points[4].mass = 4f; 

            new_simple_rig.AddBone("arm_r", 0, 1);
            new_simple_rig.bones[new_simple_rig.bones.Count-1].length[1] = measured_arm_length;
            new_simple_rig.bones[new_simple_rig.bones.Count-1].length[0] *= 0.4f; // Allow arm to flex
            new_simple_rig.AddBone("arm_l", 2, 3);
            new_simple_rig.bones[new_simple_rig.bones.Count-1].length[1] = measured_arm_length;
            new_simple_rig.bones[new_simple_rig.bones.Count-1].length[0] *= 0.4f;
            new_simple_rig.AddBone("tri_top", 0, 2);
            new_simple_rig.AddBone("tri_r", 0, 4);
            new_simple_rig.AddBone("tri_l", 2, 4);
        }
        
        // Set up full-body IK particles and bones
        complete.AddPoint(shoulder.position, "shoulder_r");
        complete.AddPoint(grip.position, "hand_r");
        complete.AddPoint((shoulder.position+Vector3.right * (neck.position[0] - shoulder.position[0])*2f), "shoulder_l");
        complete.AddPoint((grip.position+Vector3.right * (neck.position[0] - grip.position[0])*2f), "hand_l");
        complete.AddPoint(new float3(neck.position[0], hip.position[1], neck.position[2]), "body");
        complete.AddPoint(head.position, "head");
        complete.AddPoint(neck.position, "neck");
        complete.AddPoint(stomach.position, "stomach"); // 7
        complete.AddPoint(pelvis.position, "hip"); // 8
        complete.AddPoint(groin.position, "groin");
        complete.AddPoint(hip.position, "hip_r");
        complete.AddPoint(foot.position, "foot_r");
        complete.AddPoint(hip.position+Vector3.right * (neck.position[0] - hip.position[0])*2f, "hip_l");
        complete.AddPoint(foot.position+Vector3.right * (neck.position[0] - foot.position[0])*2f, "foot_l");
        
        complete.AddBone("arm_r", 0, 1);
        complete.bones[complete.bones.Count-1].length[1] = measured_arm_length;
        complete.bones[complete.bones.Count-1].length[0] *= 0.4f;
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
        
        // Create random branch 'terrain'
        int num_segments = 40;
        float x = 0;
        float y = 0;
        for(int i=0; i<num_segments+1; ++i){
            branches.AddPoint(new float3(x,y,0), "branch");
            x += UnityEngine.Random.Range(2.0f, 6.0f);
            //y += UnityEngine.Random.Range(-3.0f, 3.0f);
            //y = math.clamp(y, -2.5f, 2.5f); // Make sure we stay on screen
        }
        for(int i=0; i<num_segments; ++i){
            branches.AddBone("branch", i, i+1);
        }
        
        // Delete visible points so we don't see it when playing game
        Destroy(root.gameObject);
    }
        
    // Use law of cosines to find angles of triangle
    static float GetAngleGivenSides(float a, float b, float c){
        var top = (c*c - a*a - b*b);
        var divisor = (-2*a*b);
        if(divisor == 0f){
            return 0f;
        }
        return math.acos(math.clamp(top / divisor, -1f, 1f));
    }

    // Solve two bone IK problems
    // Arguments when dealing with arms:
    //     start_id = shoulder
    //     end_id   = hand
    //     forward actually points backwards
    //     arm_ik   = original shoulder, elbow and grip positions
    //     top      = visual bicep
    //     bottom   = visual forearm
    //     old_axis = elbow's bind space axis of rotation
    //     axis     = elbow's world space axis of rotation
    static void ApplyTwoBoneIK(int start_id, 
                               int end_id, 
                               float3 forward, 
                               float3[] ik, 
                               DisplayBone top, 
                               DisplayBone bottom, 
                               List<Verlet.Point> points, 
                               float3 old_axis, 
                               float3 axis)
    {
        var start = points[start_id]; // shoulder
        var end = points[end_id];     // hand
            
        // Get sides of triangle formed by upper and lower limb
        float dist_a =     math.distance(ik[0], ik[1]);       // Distance from bind shoulder to bind elbow
        float dist_b =     math.distance(ik[1], ik[2]);       // Distance from bind elbow to bind grip
        float dist_c =     math.distance(start.pos, end.pos); // Distance from complete.shoulder to complete.hand
        float old_dist_c = math.distance(ik[0], ik[2]);       // Distance from bind shoulder to bind hand

        // Note that dist_a and dist_b are always the same because of the length constraint between shoulder-elbow and elbow-hand
        // The only distance that can change is the shoulder-hand distance
        // That's why the computations below work

        // Get angles of triangle
        /*
                           /* <- shoulder
                          /Y|
             dist_a ->   /  |
                        /   |
             elbow  -> *X   | <- dist_c
                        \   |
             dist_b ->   \  |
                          \ |
                           \* <- hand

             X = hinge angle
             Y = base angle
        */
        var old_hinge_angle = GetAngleGivenSides(dist_a,     dist_b, old_dist_c); // Bind angle between bicep and forearm
        var hinge_angle     = GetAngleGivenSides(dist_a,     dist_b, dist_c);     // Angle between bicep and forearm
        var old_base_angle  = GetAngleGivenSides(old_dist_c, dist_a, dist_b);     // Bind angle between bicep and shoulder-hand line
        var base_angle      = GetAngleGivenSides(dist_c,     dist_a, dist_b);     // Angle between bicep and shoulder-hand line

        // Apply rotation of entire arm (shoulder->hand)
        // With only this rotation the arms are stiff (i.e. they rotate around the shoulder, but they don't bend at the elbows)
        // (end.pos - start.pos) = (hand - shoulder) = points downwards
        // forward points backwards
        var base_rotation = Quaternion.LookRotation(end.pos - start.pos, forward) * Quaternion.Inverse(Quaternion.LookRotation(end.bind_pos - start.bind_pos, Vector3.forward));
        // Apply additional rotation from IK
        base_rotation = Quaternion.AngleAxis(base_angle * Mathf.Rad2Deg, axis) * base_rotation * Quaternion.Inverse(Quaternion.AngleAxis(old_base_angle * Mathf.Rad2Deg, old_axis));

        // Apply base and hinge rotations to actual display bones
        // bind_bicep + (shoulder - bind_shoulder) = move bicep into position below the shoulder
        top.transform.position = top.bind_pos + (start.pos - start.bind_pos);
        // orient the bicep
        top.transform.rotation = base_rotation * top.bind_rot;

        // Move forearm into position below the arm
        bottom.transform.position = top.transform.position + top.transform.rotation * Quaternion.Inverse(top.bind_rot) * (bottom.bind_pos - top.bind_pos);
        bottom.transform.rotation = Quaternion.AngleAxis(hinge_angle * Mathf.Rad2Deg, axis) * base_rotation * Quaternion.Inverse(Quaternion.AngleAxis(old_hinge_angle * Mathf.Rad2Deg, old_axis)) * bottom.bind_rot;
    }

    // Calculate bone transform that matches orientation of top and bottom points, and looks in the character "forward" direction
    // start = top point
    // end   = bottom point
    void ApplyBound(DisplayBone part, float3 forward, float3 bind_forward, int start, int end){
        // Get mid point and "up" direction (from start to end point)
        // The "up" direction actually points down, since it goes from start (top) to end (bottom)
        // I think that's why the forward and bind_forward directions that this function receives are negated
        var up       = math.normalize(complete.points[end].pos      - complete.points[start].pos);
        var bind_up  = math.normalize(complete.points[end].bind_pos - complete.points[start].bind_pos);
        var mid      = (complete.points[end].pos      + complete.points[start].pos)      / 2.0f;
        var bind_mid = (complete.points[end].bind_pos + complete.points[start].bind_pos) / 2.0f;

        // Apply rotations
        // This quaternion has the same purpose as body_rotation
        // This code is actually the same as the code in the beginning of the Update function
        var rotation = Quaternion.LookRotation(up, forward) * Quaternion.Inverse(Quaternion.LookRotation(bind_up, bind_forward));
        part.transform.rotation = rotation * part.bind_rot;
        part.transform.position = mid + (float3)(rotation * (part.bind_pos - bind_mid));
    }

    // Get height of branch at given x coordinate
    float BranchHeight(float x, int start_id, int end_id){
        var start = branches.points[start_id];
        var end = branches.points[end_id];
        float branch_t = (x-start.bind_pos[0])/(end.bind_pos[0]-start.bind_pos[0]);
        return math.lerp(start.pos[1], end.pos[1], branch_t);
    }
    
    // Get height of entire branch terrain at given x coordinate
    float BranchesHeight(float x){
        for(int i=0;i<branches.bones.Count; ++i){
            var point_ids = branches.bones[i].points;
            if(x >= branches.points[point_ids[0]].pos[0] && x < branches.points[point_ids[1]].pos[0]){
                return BranchHeight(x, point_ids[0], point_ids[1]);
            }
        }
        // If not on terrain, extend horizontally forever
        if(x < 0.0f){
            return branches.points[0].pos[1];
        } else {
            return branches.points[branches.points.Count-1].pos[1];
        }
    }

    static void DrawSystem(MovementSystem system, Color color) {
        system.simple_rig.DrawBones(color); 
        for(int i=2; i<4; ++i){
            DebugDraw.Sphere(system.limb_targets[i], color, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
        }
    }

    // Prepare to draw next frame
    void Update() {
        // At this point the step function has already been called, so the display rig has been updated
        // The point of the Update function is to update the complete rig and then the display_body

        // Here we update the complete rig using the simple rig
        { // Use "arms" rig to drive full body IK rig
            var points = display.simple_rig.points;

            // Calculate midpoint and orientation of body triangle
            // Center of chest's triangle
            var bind_mid     = (points[0].bind_pos + points[2].bind_pos + points[4].bind_pos) / 3.0f;
            var mid          = (points[0].pos     +  points[2].pos      + points[4].pos)      / 3.0f;
            // Forward vector of chest's triangle
            var forward      = math.normalize(math.cross(points[0].pos      - points[2].pos,      points[0].pos      - points[4].pos));
            var bind_forward = math.normalize(math.cross(points[0].bind_pos - points[2].bind_pos, points[0].bind_pos - points[4].bind_pos));
            // Up vector that points from body to mid point between the shoulders
            var up           = math.normalize((points[0].pos      + points[2].pos)      / 2.0f - points[4].pos);
            var bind_up      = math.normalize((points[0].bind_pos + points[2].bind_pos) / 2.0f - points[4].bind_pos);
        
            // Copy hand and shoulder positions from simple rig
            for(int i=0; i<4; ++i){
                complete.points[i].pos = points[i].pos;
                complete.points[i].pinned = true;
            }

            // Starting with a simple rig (a triangle with two sticks for the arms),
            // we want to create a complete rig (a spine with two sticks for the arms and two sticks for the legs)
            // The strategy is the following: move the mid point of the chest's triangle along (bind_pos[i] - bind_mid) for every point in the complete rig
            // and rotate those vectors to match the orientation of the chest's triangle

            // body_rotation is the quaternion that rotates (bind_pos[i] - bind_mid) so that it has the orientation of the chest's triangle
            // I see it as transforming the vector to the inverse bind space, and then back to world space with the chest's orientation
            // lookRot(fwd, up) * lookRot(bindFwd, bindUp)^-1 * bindVec
            var body_rotation = math.mul(quaternion.LookRotation(forward, up), math.inverse(quaternion.LookRotation(bind_forward, bind_up)));

            // Set up spine, head and leg positions based on body rotation
            // 5  = head
            // 6  = neck
            // 7  = stomach
            // 8  = pelvis (hip)
            // 9  = groin
            // 10 = hip right
            // 11 = foot right
            // 12 = hip left
            // 13 = foot left
            for(int i=5; i<14; ++i){
                complete.points[i].pos = mid + math.mul(body_rotation, (complete.points[i].bind_pos - bind_mid));
                complete.points[i].pinned = true;
            }
            
            // Apply body compression
            // 7 = stomach
            // 8 = pelvis (hip)
            // 9 = groin
            complete.points[7].pinned = false;
            complete.points[8].pinned = false;
            var old_hip = complete.points[9].pos; // This should be called old_groin
            for(int i=7; i<=9; ++i){
                // Lerp between points 7, 8 and 9 and point 6, which is the neck
                // The higher the body compression, the closer points 7, 8 and 9 are to the neck
                complete.points[i].pos = math.lerp(complete.points[i].pos, complete.points[6].pos, body_compress_amount);
            }
            // Move the stomach and the pelvis (hip) back
            complete.points[7].pos -= forward * body_compress_amount * 0.2f;
            complete.points[8].pos -= forward * body_compress_amount * 0.2f;

            // Move the legs towards the new position of the groin
            for(int i=10; i<14; ++i){
                complete.points[i].pos += complete.points[9].pos - old_hip;
            }

            // Move feet to foot targets
            for(int i=0; i<2; ++i){
                complete.points[11+i*2].pos = display.limb_targets[2+i];
            }
            
            // Enforce bone length constraints
            for(int i=0; i<2; ++i){
                complete.EnforceDistanceConstraints();
            }
        }

        // Here we update the display_body using the complete rig
        { // Apply full body IK rig to visual deformation bones
            var points = complete.points;

            // Get torso orientation and position
            // Note that here we are working with the complete rig, not the simple rig
            // 0 = right shoulder
            // 2 = left shoulder
            // 9 = groin
            // Center of torso
            var bind_mid     = (points[0].bind_pos + points[2].bind_pos + points[9].bind_pos) / 3.0f;
            var mid          = (points[0].pos      + points[2].pos      + points[9].pos)      / 3.0f;
            // Forward vector of torso negated, so it points backwards
            var forward      = -math.normalize(math.cross(points[0].pos      - points[2].pos,      points[0].pos -      points[9].pos));
            var bind_forward = -math.normalize(math.cross(points[0].bind_pos - points[2].bind_pos, points[0].bind_pos - points[9].bind_pos));
            // Up vector that points from groin to mid point between the shoulders
            var up           = math.normalize((points[0].pos      + points[2].pos)/2.0f      - points[9].pos);
            var bind_up      = math.normalize((points[0].bind_pos + points[2].bind_pos)/2.0f - points[9].bind_pos);
        
            // Apply core bones
            // 5 = head
            // 6 = neck
            ApplyBound(display_body.head, forward, bind_forward, 5, 6);
            // 6 = neck
            // 7 = stomach
            ApplyBound(display_body.chest, forward, bind_forward, 6, 7);
            // 7 = stomach
            // 8 = pelvis (hip)
            ApplyBound(display_body.belly, forward, bind_forward, 7, 8);
            // 8 = pelvis (hip)
            // 9 = groin
            ApplyBound(display_body.pelvis, forward, bind_forward, 8, 9);

            // Arm IK
            for(int i=0; i<2; ++i)
            {
                var top = display_body.arm_top_r; // Bicep
                var bottom = display_body.arm_bottom_r; // Forearm
                if (i==1)
                {
                    top = display_body.arm_top_l; // Bicep
                    bottom = display_body.arm_bottom_l; // Forearm
                }

                // 0 = right shoulder
                // 1 = right hand
                // 2 = left shoulder
                // 3 = left hand
                int start_id = i * 2;            // Shoulder
                int end_id   = i * 2 + 1;        // Hand
                var start    = points[start_id]; // Shoulder point
                var end      = points[end_id];   // Hand point

                // Adjust elbow target position
                float ik_driver = 1.0f;
                var ik_forward_amount = -ik_driver * 0.8f;
                var ik_up_amount = 0.1f + ik_driver * 0.5f;
                // Start at mid point between the shoulders, move up by ik_up_amount and move forward by ik_forward_amount
                var elbow_point      = ((points[2].pos      + points[0].pos)      * 0.5f + up      * ik_up_amount + forward      * ik_forward_amount);
                var bind_elbow_point = ((points[2].bind_pos + points[0].bind_pos) * 0.5f + bind_up * ik_up_amount + bind_forward * ik_forward_amount);

                if(debug_info.draw_elbow_ik_target)
                {
                    DebugDraw.Line((start.pos + end.pos) * 0.5f, elbow_point, Color.red, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                    DebugDraw.Sphere(elbow_point, Color.red, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                }

                // ((mid point between hand and shoulder) - elbow IK target) cross (shoulder - hand)
                // This looks to me like the axis of rotation of the elbow
                // It points outwards from the body (remember to use your left hand when computing cross products in a left-handed coordinate system like Unity's)
                // old_axis is bind space
                // axis is world space
                var old_axis = math.normalize(math.cross((end.bind_pos + start.bind_pos) * 0.5f - bind_elbow_point, start.bind_pos - end.bind_pos));
                var axis     = math.normalize(math.cross((end.pos      + start.pos)      * 0.5f - elbow_point,      start.pos      - end.pos));

                // start_id = shoulder
                // end_id   = hand
                // forward actually points backwards
                // arm_ik   = original shoulder, elbow and grip positions
                // top      = visual bicep
                // bottom   = visual forearm
                // old_axis = elbow's bind space axis of rotation
                // axis     = elbow's world space axis of rotation
                ApplyTwoBoneIK(start_id, end_id, forward, arm_ik, top, bottom, complete.points, old_axis, axis);
                
                if(debug_info.draw_ik_final)
                {
                    DebugDraw.Line(points[start_id].pos, bottom.transform.position, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                    DebugDraw.Line(points[end_id].pos, bottom.transform.position, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                }
            }

            // Leg IK
            for(int i=0; i<2; ++i){
                var top = display_body.leg_top_r; // Thigh
                var bottom = display_body.leg_bottom_r; // Calf
                if(i==1){
                    top = display_body.leg_top_l; // Thigh
                    bottom = display_body.leg_bottom_l; // Calf
                }

                // 10 = right hip
                // 11 = right foot
                // 12 = left hip
                // 13 = left foot
                int start = i*2+10; // Hip
                int end = i*2+1+10; // Foot

                // Points downwards from hip to foot
                var leg_dir = points[end].pos - points[start].pos;

                // Get knee direction
                var leg_dir_flat = math.normalize(new float2(math.dot(leg_dir, forward), math.dot(leg_dir, up)));
                var leg_forward = leg_dir_flat[0] * up + leg_dir_flat[1] * -forward;
                //DebugDraw.Line(mid, mid + leg_forward, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);

                // Get base whole-leg rotation
                var bind_rotation = Quaternion.LookRotation(points[end].bind_pos - points[start].bind_pos, Vector3.forward);
                var rotation = Quaternion.LookRotation(leg_dir, leg_forward) * bind_rotation;

                // Get knee bend axis
                var old_axis = bind_rotation * Vector3.right;
                var axis = rotation * Vector3.right;

                ApplyTwoBoneIK(start, end, leg_forward, leg_ik, top, bottom, complete.points, old_axis, axis);
                
                if(debug_info.draw_ik_final){
                    DebugDraw.Line(points[start].pos, bottom.transform.position, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                    DebugDraw.Line(points[end].pos, bottom.transform.position, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                }
            }

            // Head look
            // head_look_y: 50 = max look down, -70 = max look up
            // head_look_x: -90 to 90

            // Get head target in head transform space
            var target = math.normalize(display_body.head.transform.InverseTransformPoint(look_target));
            // Using sin here is not correct (should be asin or something), but looks ok so keeping it for now
            var head_look_y = math.sin(target.x)*Mathf.Rad2Deg;
            // Flatten look direction to solve other rotation axis
            var temp = target;
            temp.x = 0.0f;
            temp = math.normalize(temp);
            var head_look_x = -math.sin(temp.y)*Mathf.Rad2Deg;

            // Apply head transform
            display_body.head.transform.rotation = display_body.head.transform.rotation * Quaternion.AngleAxis(head_look_x, new Vector3(1.0f, 0.0f, 0.0f)) * Quaternion.AngleAxis(head_look_y, new Vector3(0.0f, 1.0f, 0.0f));
            if(head_look_y > 0.0f){
                // It's important not to take the scale into account here, which is why we use TransformDirection instead of TransformVector
                display_body.head.transform.position = display_body.head.transform.position + (Vector3)((display_body.head.transform.TransformDirection(new Vector3(1.0f, 0.0f, 0.0f))) * head_look_y * -0.001f);
            }
            
            if(debug_info.draw_head_look){
                DebugDraw.Sphere(look_target, Color.red, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
                DebugDraw.Line(display_body.head.transform.position, look_target, Color.red, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
            }
        }

        // Debug draw skeleton
        branches.DrawBones(new Color(0.5f, 0.5f, 0.1f, 1.0f));
        if(debug_info.draw_walk_rig){ DrawSystem(walk, Color.red); }
        if(debug_info.draw_display_simple_rig){ DrawSystem(display, Color.white); }
        if(debug_info.draw_display_complete_rig){ complete.DrawBones(Color.white); }
        if(debug_info.draw_ik_final){
            for(int i=2; i<complete.bones.Count-2; ++i){
                DebugDraw.Line(complete.points[complete.bones[i].points[0]].pos,complete.points[complete.bones[i].points[1]].pos, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
            }       
            DebugDraw.Line(complete.points[0].pos,complete.points[2].pos, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
            DebugDraw.Line(complete.points[10].pos,complete.points[12].pos, Color.white, DebugDraw.Lifetime.OneFrame, DebugDraw.Type.Xray);
        }
        display_gibbon.SetActive(debug_info.draw_gibbon);

        debug_info.DrawWindow();

        if(Input.GetKeyDown(KeyCode.Tab)){
            Time.timeScale = (Time.timeScale == 1.0f)?0.1f:1.0f;
        }
    }
    
    float MoveTowardsF(float a, float b, float max_dist){
        float len = math.distance(a,b);
        if(len < max_dist){
            return b;
        } else {
            return a + (b-a)/len*max_dist;
        }
    }

    float3 MoveTowardsVec(float3 a, float3 b, float max_dist){
        float len = math.distance(a,b);
        if(len < max_dist){
            return b;
        } else {
            return a + (b-a)/len*max_dist;
        }
    }
    
    void Swap(ref float3 a, ref float3 b){
        var temp = a;
        a = b;
        b = temp;
    }
    
    void PreventHandsFromCrossingBody(Verlet.System rig) {    
        for(int i=0; i<2; ++i){
            // For the right hand this vector points from the left shoulder to the right shoulder
            // For the left hand it points in the opposite direction
            var side_dir = math.normalize(rig.points[0].pos - rig.points[2].pos) * (1-i*2);
            float shoulder_d = math.dot(rig.points[i*2].pos, side_dir);
            float hand_d = math.dot(rig.points[i*2+1].pos, side_dir);
            float new_d = math.max(hand_d, shoulder_d);
            //if (i == 0) { DebugDraw.Line(new float3(0.0f, 0.0f, 0.0f), side_dir, Color.red, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray); }
            //if (i == 0) { DebugDraw.Line(new float3(0.0f, 0.0f, 0.0f), rig.points[i * 2].pos, Color.green, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray); }
            //if (i == 0) { DebugDraw.Line(new float3(0.0f, 0.0f, 0.0f), rig.points[i * 2 + 1].pos, Color.blue, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray); }

            //if (i == 0) { DebugDraw.Sphere(rig.points[i * 2 + 1].pos, Color.red, Vector3.one * 0.05f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray); }
            //if (i == 0) { DebugDraw.Line(rig.points[i * 2 + 1].pos, rig.points[i * 2 + 1].pos + (new_d - hand_d) * side_dir * 5.0f, Color.red, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray); }
            rig.points[i*2+1].pos += (new_d - hand_d) * side_dir;
            //if (i == 0) { DebugDraw.Sphere(rig.points[i * 2 + 1].pos, Color.green, Vector3.one * 0.05f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray); }
        }
    }

    // Apply actual controls and physics
    void Step(float step) {
        // Transform controls to axes
        float horz_input = 0f;
        if(Input.GetKey(KeyCode.D)){
            horz_input = 1f; // Move to the right (along the +X direction)
        }
        if(Input.GetKey(KeyCode.A)){
            horz_input = -1f; // Move to the left (along the -X direction)
        }

        // Max speed of 7 m/s while running
        float max_speed = 7.0f;
        
        // Player can influence horizontal velocity
        simple_vel[0] += horz_input * step * 5f;
        simple_vel[0] = math.clamp(simple_vel[0], -max_speed, max_speed);
        
        // Don't allow speed < 1.0 m/s, don't need to worry about idle animations in an endless runner
        if(horz_input == 0f && math.abs(simple_vel[0]) < 1.0f){
            simple_vel[0] = MoveTowardsF(simple_vel[0], simple_vel[0] >= 0.0f ? 1.0f : -1.0f, step);
        }

        // Smooth out vertical position on branch by checking height forwards and back
        // Smoothing helps make the motion look good at the points where the slope of the terrain changes abruptly
        // It makes the height of the character not match the height of the branches at those points
        // In points like this one: /\ -> The height of the character is lower than the corner
        // In points like this one: \/ -> the height of the character is higher than the corner
        var future_pos = simple_pos + simple_vel * 0.1f;
        future_pos[1] = BranchesHeight(future_pos[0]);
        var past_pos = simple_pos + simple_vel * -0.1f;
        past_pos[1] = BranchesHeight(past_pos[0]);
        var smoothed_pos = (future_pos + past_pos + simple_pos) / 3.0f;

        // Get slope and use it to modify running speed
        var slope_vec = math.normalizesafe(future_pos - simple_pos);
        float slope_speed_mult = math.abs(slope_vec[0]);

        // Apply modified running speed to position
        var effective_vel = simple_vel * slope_speed_mult;
        simple_pos += effective_vel * step;

        simple_pos[1] = BranchesHeight(simple_pos[0]);
        simple_vel[1] = 0.0f;

        // If on ground, look in the direction you are moving
        // Here we calculate the forward vector of the chest's triangle
        var forward = math.normalize(math.cross(display.simple_rig.points[0].pos - display.simple_rig.points[2].pos, display.simple_rig.points[0].pos - display.simple_rig.points[4].pos));
        look_target = (float3)display_body.head.transform.position + forward * 0.1f;
        look_target += future_pos - past_pos;

        if(debug_info.draw_smoothing){
            DebugDraw.Sphere(future_pos, Color.red, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            DebugDraw.Sphere(past_pos, Color.blue, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            DebugDraw.Sphere(smoothed_pos, Color.green, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            DebugDraw.Sphere(simple_pos, Color.yellow, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            DebugDraw.Line(future_pos, smoothed_pos, Color.green, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            DebugDraw.Line(past_pos, smoothed_pos, Color.green, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
        }

        { // Run animation
            // Vary between different gaits based on speed and time
            // y\ =\sin(x\cdot2.3)+\sin(x\cdot1.7)\ \left\{0<y\ <\ 1\right\}\ \left\{x\ >0\right\}
            quad_amount = math.clamp((math.sin(Time.time * 2.3f) + math.sin(Time.time * 1.7f)), 0.0f, 1.0f);

            // Determine how far to lean forwards
            // y\ =\ \sin(x)*0.2+0.3\ \left\{x\ >0\right\}
            var walk_lean = math.sin(Time.time) * 0.2f + 0.3f;
            var lean = walk_lean;

            // Adjust stride frequency based on speed
            // y\ =\ 8\ /\ (\pi\ *\ 2)\ *\ (\operatorname{abs}(x)+1.0)^{0.4}\ \left\{1\ <\ x\ <\ 7\right\}
            float speed_mult = 8f / (math.PI * 2f) * math.pow((math.abs(effective_vel[0]) + 1.0f), 0.4f);
            // When the character moves faster, we increase the walk_time faster
            walk_time += step*speed_mult;

            // Compress body during quadruped gallop
            walk.body_compress_amount = 0.0f;

            // Adjust COM height based on gait
            var target_com = simple_pos;
            target_com[1] = smoothed_pos[1];
            // walk_height is what raises the COM from the floor. Notice how it draws a sine wave
            var walk_height = base_walk_height + math.sin((walk_time+0.25f) * math.PI * 4.0f) * math.abs(effective_vel[0]) * 0.015f / speed_mult + math.abs(effective_vel[0]) * 0.01f;
            target_com[1] += walk_height;
            // The target COM moves up and down. The lower it is (i.e. the closer it is to simple_pos), the more the character leans forward
            //DebugDraw.Sphere(target_com, Color.red, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            //DebugDraw.Sphere(simple_pos, Color.green, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            target_com[1] = math.lerp(target_com[1], simple_pos[1], math.abs(lean)*0.15f);
            //DebugDraw.Sphere(target_com, Color.blue, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);

            if (debug_info.draw_smoothing)
            {
                DebugDraw.Sphere(target_com, Color.cyan, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            }

            // Get ground slope again for use later
            var left = simple_pos - new float3(0.1f,0.0f,0.0f);
            var right = simple_pos + new float3(0.1f,0.0f,0.0f);
            left[1] = BranchesHeight(left[0]);
            right[1] = BranchesHeight(right[0]);
            // move_dir always points to the right (i.e. down the +X axis)
            float3 move_dir = math.normalize(right - left);
            //DebugDraw.Line(simple_pos + new float3(0.0f, 0.5f, 0.0f), (simple_pos + new float3(0.0f, 0.5f, 0.0f)) + (move_dir), Color.red, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
            //DebugDraw.Sphere((simple_pos + new float3(0.0f, 0.5f, 0.0f)) + (move_dir), Color.red, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);

            { // Simulate the walk simple rig
                var rig = walk.simple_rig;

                // Apply gravity and velocity to all the points in the walk simple rig
                rig.StartSim(step);

                for(int j=0; j<4; ++j){
                    // Calculate the COM
                    float total_mass = 0f;
                    var com = float3.zero;
                    for(int i=0; i<rig.points.Count; ++i){
                        // 1 is the right hand and 3 is the left hand
                        // We don't want to consider those points when calculating the COM
                        // We only consider the points that make up the chest's triangle (0 - right should, 2 - left shoulder and 4 - body)
                        if(i!=1 && i!=3){
                            com += rig.points[i].pos * rig.points[i].mass;
                            total_mass += rig.points[i].mass;
                        }
                    }

                    // Adjust all free points to match target COM
                    com /= total_mass;
                    var offset = (float3)target_com - com;
                    for(int i=0; i<rig.points.Count; ++i){
                        if(i!=1 && i!=3){
                            rig.points[i].pos += offset * 0.2f;
                        }
                    }

                    // Apply torque to keep torso upright and forward-facing
                    float step_sqrd = step*step;
                    float force = 20f;
                    // Calculate the forward vector of the chest's triangle
                    var forward2 = math.normalize(math.cross(rig.points[0].pos - rig.points[2].pos, rig.points[0].pos - rig.points[4].pos));
                    var flat_forward = math.normalize(new float3(forward2[0],0,forward2[2]));
                    // flat_foward is simply the forward vector of the chest's triangle with no Y component, so it's always perfectly aligned with the +X and -X axes
                    // Multiplying flat_forward by lean causes the length of the vector to shrink when the character is upright, and to grow when the character is leaning
                    // Adding (0, 1, 0) just makes the vector point upwards diagonally 
                    float3 top_force = (lean*flat_forward + new float3(0,1,0)) * force;
                    //DebugDraw.Line(simple_pos + new float3(0.0f, 0.5f, 0.0f), (simple_pos + new float3(0.0f, 0.5f, 0.0f)) + (flat_forward), Color.red, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
                    //DebugDraw.Line(simple_pos + new float3(0.0f, 0.5f, 0.0f), (simple_pos + new float3(0.0f, 0.5f, 0.0f)) + (lean * flat_forward), Color.red, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
                    //DebugDraw.Line(simple_pos + new float3(0.0f, 0.5f, 0.0f), (simple_pos + new float3(0.0f, 0.5f, 0.0f)) + (lean * flat_forward + new float3(0, 1, 0)), Color.red, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
                    //DebugDraw.Line(simple_pos + new float3(0.0f, 0.5f, 0.0f), (simple_pos + new float3(0.0f, 0.5f, 0.0f)) + (top_force), Color.red, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
                    // These three lines are what keep the chest upright
                    // The torque is carefully balanced here
                    rig.points[4].pos += step_sqrd * -top_force;                 // Body
                    rig.points[0].pos += step_sqrd * top_force * 0.5f;           // Right shoulder
                    rig.points[2].pos += step_sqrd * top_force * 0.5f;           // Left shoulder
                    // These two lines are what make the character look forward when it moves. Without them, the character walks sideways
                    rig.points[0].pos[2] -= step_sqrd * effective_vel[0] * 2.0f; // Right shoulder
                    rig.points[2].pos[2] += step_sqrd * effective_vel[0] * 2.0f; // Left shoulder
                    
                    // Add rotational force to body if needed
                    // This is what makes the chest sway from side to side (i.e. rotate around the Z axis)
                    for(int i=0; i<2; ++i){
                        var walk_rotate = (math.cos((walk_time + tilt_offset) * math.PI * 2.0f + math.PI * i)) * 0.2f;
                        var rotate = walk_rotate;
                        // Points 0 and 2 are modified here (right shoulder and left shoulder, respectively)
                        rig.points[i*2].pos[0] += step_sqrd * -3.0f * rotate * effective_vel[0] / speed_mult;
                    }
                    
                    // Move arms out to sides
                    float speed = math.abs(effective_vel[0]) / max_speed;
                    for(int i=0; i<2; ++i){
                        // This is what causes the arms to flail
                        var arms_up = math.abs(speed * (math.sin(Time.time * ((i==1)?2.5f:2.3f))*0.3f+0.7f));
                        // Points 1 and 3 are modified here (right hand and left hand, respectively)
                        // They are moved along the vector that connects the shoulders (points 0 and 2)
                        rig.points[1+i*2].pos += step_sqrd * (rig.points[0].pos - rig.points[2].pos) * (1.5f + speed * 2.0f + arms_up * 2.0f) * (1 - i * 2) * 2f;
                        rig.points[1+i*2].pos[1] += step_sqrd * 10.0f * arms_up  * arms_up;
                        // Bone 0 is the right arm, while bone 1 is the left arm
                        // Here we update the maximum length
                        rig.bones[i].length[1] = rig.bones[0].length[0] / 0.4f * (math.lerp(0.95f, 0.8f, math.min(speed*0.25f, 1.0f) + math.sin(arms_up * math.PI)*0.1f));
                    }
                    
                    // This prevents the arms from colliding with the body and penetrating it when the character turns
                    // TODO: I still don't understand the math behind this function
                    PreventHandsFromCrossingBody(rig);
                    
                    // Make sure hands don't go through floor
                    for(int i=0; i<2; ++i){
                        // Points 1 and 3 are modified here (right hand and left hand, respectively)
                        rig.points[i*2+1].pos[1] = math.max(rig.points[i*2+1].pos[1], BranchesHeight(rig.points[i*2+1].pos[0]));
                    }
                    
                    for(int i=0; i<2; ++i){
                        rig.EnforceDistanceConstraints();
                    }
                }

                rig.EndSim();
                
                // Calculate leg targets
                // 2 is right foot
                // 3 is left foot
                for(int i=0; i<2; ++i){
                    var offset = math.lerp(gallop_offset, quad_gallop_offset, quad_amount);
                    walk.limb_targets[2+i] = simple_pos;
                    // This vector points in the movement direction. Note how it's length depends on the effective velocity. That causes the strides to shorten or lengthen depending on the velocity of the character
                    walk.limb_targets[2+i] += (move_dir * (math.cos(walk_time * math.PI * 2.0f + math.PI*i))*0.2f) * effective_vel[0] / speed_mult;
                    // This vector points from the left shoulder to the right shoulder for the right foot, and in the opposite direction for the left foot
                    walk.limb_targets[2+i] += (rig.points[0].pos - rig.points[2].pos) * (1.0f-2.0f*i) * (0.3f);
                    // The height of the target matches the height of the branch
                    // It float aboves the ground in front and behind the character, and it sinks in the middle
                    walk.limb_targets[2+i][1] = BranchesHeight(walk.limb_targets[2+i][0]);
                    walk.limb_targets[2+i][1] += (-math.sin(walk_time * math.PI * 2.0f + math.PI*i) + 1.0f)*0.2f * (math.pow(math.abs(effective_vel[0]) + 1.0f, 0.3f) - 1.0f) * (1.0f);

                    //DebugDraw.Sphere(walk.limb_targets[2 + i], (i == 0) ? Color.green : Color.red, Vector3.one * 0.05f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
                }
            }
        } 
        
        { // Combine source rigs into display rig
            // Start calculating COM line if needed
            var old_com = float3.zero;
            if(debug_info.draw_com_line) {
                float total_mass = 0.0f;
                var points = display.simple_rig.points;
                for(int i=0; i<points.Count; ++i){
                    if(points[i].pinned == false){
                        old_com += points[i].pos * points[i].mass;
                        total_mass += points[i].mass;
                    }
                }
                old_com /= total_mass;
            }

            // Interpolate between source rigs
            // Since we only support running right now, there is no interpolation
            // The walk right is simply copied into the display rig
            for(int i=0; i<display.simple_rig.points.Count; ++i){
                display.simple_rig.points[i].old_pos = display.simple_rig.points[i].pos;
                display.simple_rig.points[i].pos = walk.simple_rig.points[i].pos;
            }
            for(int i=0; i<2; ++i){
                display.simple_rig.EnforceDistanceConstraints();
            }
            for(int i=0; i<4; ++i){
                display.limb_targets[i] = walk.limb_targets[i];
            }
            body_compress_amount = walk.body_compress_amount;
            
            // Draw COM line
            if(debug_info.draw_com_line) {
                var com = float3.zero;
                {
                    float total_mass = 0.0f;
                    var points = display.simple_rig.points;
                    for(int i=0; i<points.Count; ++i){
                        if(points[i].pinned == false){
                            com += points[i].pos * points[i].mass;
                            total_mass += points[i].mass;
                        }
                    }
                    com /= total_mass;
                }
                debug_info.com_lines.Add(DebugDraw.Line(old_com, com, Color.green, DebugDraw.Lifetime.Persistent, DebugDraw.Type.Xray));
            }
        }

        if(debug_info.draw_simple_point){
            DebugDraw.Sphere(simple_pos, Color.yellow, Vector3.one * 0.1f, Quaternion.identity, DebugDraw.Lifetime.OneFixedUpdate, DebugDraw.Type.Xray);
        }

        // Move game camera to track character
        {
            var cam_pos = Camera.main.transform.position;
            // Get COM
            float total_mass = 0.0f;
            var com = new float3(0.0f, 0.0f, 0.0f);
            var points = display.simple_rig.points;
            for(int i=0; i<points.Count; ++i){
                com += points[i].pos * points[i].mass;
                total_mass += points[i].mass;
            }
            com /= total_mass;
            // Track COM position
            cam_pos[0] = com[0] + simple_vel[0] * 0.1f;
            cam_pos[1] = com[1];
            Camera.main.transform.position = cam_pos;
        }

    }

    private void FixedUpdate() {
        Step(Time.fixedDeltaTime);
    }
}
}