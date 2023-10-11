#version 450
#extension GL_EXT_nonuniform_qualifier : enable

layout(set = 0, binding = 0) uniform UniformBufferCamera 
{
	vec3 camPos;
    mat4 view;
    mat4 proj;
} ubc;

layout(set = 1, binding = 0) uniform sampler2D texSampler[];

layout(set = 2, binding = 0) uniform UniformBufferModel
{
    mat4 model;
} ubo[];


layout(set = 3, binding = 0) uniform UniformBufferMaterial
{
	vec3  albedo;
	vec2 offset;
	vec2 tilling;
	float metallic;
    float roughness;	
	float normal;
	float ao;
	uint albedoMap;
	uint normalMap;
	uint metallicMap;
    uint roughnessMap;
	uint aoMap;
	uint castShadow;
} ubm[];

layout(set = 4, binding = 0) uniform UniformBufferLight
{
	vec3 position;
    vec3 color;
	vec3 direction;
	float range;
	float spotAngle;
	uint status;//DirLight = 0 ; PointLight = 1 ; SpotLight = 2
	uint shadowID;
} ubl[];

layout(set = 5, binding = 0) uniform UniformBufferDiver
{
	uint maxLight;
	uint maxShadow;
	float u_time;
	float gamma;
}ubd;

layout (set = 6, binding = 0) uniform samplerCube samplerCubeMap;

layout(set = 7, binding = 0) uniform sampler2D shadowSampler[];

layout(set = 8, binding = 0) uniform UniformBufferShadow
{
    mat4 view;
    mat4 proj;
    vec3 position;
} ubs[];

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;
layout(location = 2) in vec3 inColor;
layout(location = 3) in vec2 inTexCoord;
layout(location = 4) in vec3 inTangent;

//Instanced
layout(location = 5) in int index_ubo;
layout(location = 6) in int index_material;

layout(location = 0) out vec2 fragTexCoord;
layout(location = 1) out vec3 Color;
layout(location = 2) out vec3 WorldPos;
layout(location = 3) out mat3 TBN;
layout(location = 6) out int imaterial;

void main() 
{
    float amplitude = 3.5; 
    float frequency = 2.0; 
    float sizePos = 0.1;

    vec3 newPosition = inPosition + vec3(0.0, amplitude * sin(ubd.u_time * frequency + (inPosition.x + inPosition.z) * sizePos), 0.0);

    vec3 normalDisplacement = vec3(sin(ubd.u_time * frequency + (inPosition.x + inPosition.z) * sizePos) * (amplitude / 12), 0.0, sin(ubd.u_time * frequency + (inPosition.x + inPosition.z) * sizePos) * (amplitude / 12));

    // Update the normal to account for displacement
    vec3 newNormal = normalize(inNormal + normalDisplacement);
    newNormal.y = inNormal.y; // Keep the y-component unchanged

    fragTexCoord = ubm[index_material].offset + inTexCoord * ubm[index_material].tilling;
    WorldPos = vec3(ubo[index_ubo].model * vec4(newPosition, 1.0));
    Color = inColor;

    vec3 T = normalize(vec3(ubo[index_ubo].model * vec4(inTangent, 0.0)));
    vec3 N = newNormal; 
    T = normalize(T - dot(T, N) * N);
    vec3 B = cross(N, T);
    TBN = mat3(T, B, N);

    imaterial = index_material;
    gl_Position = ubc.proj * ubc.view * ubo[index_ubo].model * vec4(newPosition, 1.0);
}