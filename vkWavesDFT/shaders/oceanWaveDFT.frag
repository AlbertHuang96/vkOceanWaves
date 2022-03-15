#version 450


layout (location = 0) in vec3 inHalfWay;
layout (location = 1) in float fog_factor;
layout (location = 2) in vec3 inNormal;
layout (location = 3) in vec3 inLightVec;

layout (location = 0) out vec4 outFragColor;

void main() 
{
    vec3 normal1         = normalize(inNormal);
    vec3 light_vector1   = normalize(inLightVec);
    vec3 halfway_vector1 = normalize(inHalfWay);
 
    vec4 c = vec4(1,1,1,1);
 
    vec4 emissive_color = vec4(1.0, 1.0, 1.0,  1.0);
    vec4 ambient_color  = vec4(0.0, 0.65, 0.75, 1.0);
    vec4 diffuse_color  = vec4(0.5, 0.65, 0.75, 1.0);
    vec4 specular_color = vec4(1.0, 0.25, 0.0,  1.0);
 
    float emissive_contribution = 0.00;
    float ambient_contribution  = 0.30;
    float diffuse_contribution  = 0.30;
    float specular_contribution = 1.80;
 
    float d = dot(normal1, light_vector1);
    bool facing = d > 0.0;
 
    outFragColor = emissive_color * emissive_contribution +
                   ambient_color  * ambient_contribution  * c +
                   diffuse_color  * diffuse_contribution  * c * max(d, 0) +
                    (facing ?
                   specular_color * specular_contribution * c * max(pow(dot(normal1, halfway_vector1), 120.0), 0.0) :
                   vec4(0.0, 0.0, 0.0, 0.0));
 
    outFragColor = outFragColor * (1.0 - fog_factor) + vec4(0.25, 0.75, 0.65, 1.0) * (fog_factor);
 
    outFragColor.a = 1.0;
}	
