#version 450

layout (location = 0) in vec3 inPos;
layout (location = 1) in vec3 inNormal;

layout (binding = 0) uniform UBO 
{
	mat4 projection;
	mat4 model;
	vec4 viewPos;
	vec4 lightPos;
} ubo;

layout (location = 0) out vec3 outHalfWay;
layout (location = 1) out float fog_factor;
layout (location = 2) out vec3 outNormal;
layout (location = 3) out vec3 outLightVec;

out gl_PerVertex 
{
    vec4 gl_Position;   
};

// view alrealy multiplied into model matrix
void main() 
{
	gl_Position = ubo.model * vec4(inPos.xyz, 1.0);
	fog_factor = min(-gl_Position.z / 500.0, 1.0);
	gl_Position = ubo.projection * gl_Position;

    vec4 pos = ubo.model * vec4(inPos, 1.0);
	outNormal = normalize(mat3(inverse(transpose(ubo.model))) * inNormal);

	vec3 lPos = normalize(mat3(ubo.model) * ubo.lightPos.xyz);
    outLightVec = lPos - pos.xyz;

	outHalfWay = outLightVec - lPos;
    	
}
