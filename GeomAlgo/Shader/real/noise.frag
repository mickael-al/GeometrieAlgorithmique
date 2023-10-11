#version 450
#extension GL_EXT_nonuniform_qualifier : enable

const float PI = 3.14159265359;

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
} ubd;

layout(set = 6, binding = 0) uniform samplerCube samplerCubeMap;

layout(set = 7, binding = 0) uniform sampler2D shadowSampler[];

layout(set = 8, binding = 0) uniform UniformBufferShadow
{
    mat4 view;
    mat4 proj;
    vec3 position;
} ubs[];

layout(location = 0) in vec2 fragTexCoord;
layout(location = 1) in vec3 Color;
layout(location = 2) in vec3 WorldPos;
layout(location = 3) in mat3 TBN;
layout(location = 6) in flat int imaterial;

layout(location = 0) out vec4 outColor;

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH * NdotH;

    float num = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r * r) / 8.0;

    float num = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return num / denom;
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2 = GeometrySchlickGGX(NdotV, roughness);
    float ggx1 = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

float isInShadow(vec4 fragPosLightSpace,uint id)
{
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;    
    if(projCoords.x < 0.0 || projCoords.x > 1.0 || projCoords.y < 0.0 || projCoords.y > 1.0)
    {
        return 1.0;
    }
    float shadowDepth = texture(shadowSampler[id], projCoords.xy).r;    
    shadowDepth = (shadowDepth - fragPosLightSpace.z) / (fragPosLightSpace.w - fragPosLightSpace.z);   
    if(shadowDepth == 0)
    {
        return 1.0;
    }
    return (projCoords.z-0.5) >= shadowDepth ? 0.0 : 1.0;
}

vec2 fade(vec2 t) {return t*t*t*(t*(t*6.0-15.0)+10.0);}

vec4 permute(vec4 x){return mod(((x*34.0)+1.0)*x, 289.0);}
vec4 taylorInvSqrt(vec4 r){return 1.79284291400159 - 0.85373472095314 * r;}
vec3 fade(vec3 t) {return t*t*t*(t*(t*6.0-15.0)+10.0);}

float cnoise(vec2 P,float res){
  vec4 Pi = floor(P.xyxy) + vec4(0.0, 0.0, 1.0, 1.0);
  vec4 Pf = fract(P.xyxy) - vec4(0.0, 0.0, 1.0, 1.0);
  Pi = mod(Pi, 289.0); // To avoid truncation effects in permutation
  vec4 ix = Pi.xzxz;
  vec4 iy = Pi.yyww;
  vec4 fx = Pf.xzxz;
  vec4 fy = Pf.yyww;
  vec4 i = permute(permute(ix) + iy);
  vec4 gx = 2.0 * fract(i * 0.0243902439) - 1.0; // 1/41 = 0.024...
  vec4 gy = abs(gx) - 0.5;
  vec4 tx = floor(gx + 0.5);
  gx = gx - tx;
  vec2 g00 = vec2(gx.x,gy.x);
  vec2 g10 = vec2(gx.y,gy.y);
  vec2 g01 = vec2(gx.z,gy.z);
  vec2 g11 = vec2(gx.w,gy.w);
  vec4 norm = 1.79284291400159*res - 0.85373472095314 * 
    vec4(dot(g00, g00), dot(g01, g01), dot(g10, g10), dot(g11, g11));
  g00 *= norm.x;
  g01 *= norm.y;
  g10 *= norm.z;
  g11 *= norm.w;
  float n00 = dot(g00, vec2(fx.x, fy.x));
  float n10 = dot(g10, vec2(fx.y, fy.y));
  float n01 = dot(g01, vec2(fx.z, fy.z));
  float n11 = dot(g11, vec2(fx.w, fy.w));
  vec2 fade_xy = fade(Pf.xy);
  vec2 n_x = mix(vec2(n00, n01), vec2(n10, n11), fade_xy.x);
  float n_xy = mix(n_x.x, n_x.y, fade_xy.y);
  return 2.3 * n_xy;
}

float remap(float value, float oldMin, float oldMax, float newMin, float newMax) {
    return newMin + (value - oldMin) * (newMax - newMin) / (oldMax - oldMin);
}

void main()
{   
    vec4 cr = texture(texSampler[ubm[imaterial].albedoMap], fragTexCoord).rgba;
    if(cr.a == 0.0)
    {
        discard;
        return;
    }
    float n = remap(cnoise(WorldPos.xz+vec2(ubd.u_time*0.25,0),0.4), 0.0, 1.0, 0.4, 1.0);    
    vec3 color = cr.rgb * ubm[imaterial].albedo * Color;
    vec3 ambient = vec3(0.003) * color * ubm[imaterial].ao * texture(texSampler[ubm[imaterial].aoMap], fragTexCoord).rgb;
    vec3 metallic = texture(texSampler[ubm[imaterial].metallicMap], fragTexCoord).rgb * ubm[imaterial].metallic;
    float roughness = texture(texSampler[ubm[imaterial].roughnessMap], fragTexCoord).r * ubm[imaterial].roughness;
    vec3 normal = texture(texSampler[ubm[imaterial].normalMap], fragTexCoord).rgb * 2.0;
    normal = mix(vec3(0.5, 0.5, 1.0), normal, ubm[imaterial].normal);
    normal = normalize(normal * 2.0 - 1.0);

    vec3 N = normalize(TBN * normal);
    vec3 V = normalize(ubc.camPos - WorldPos);

    vec3 F0 = vec3(0.04); 
    F0 = mix(F0, color, metallic);    

    vec3 reflectedSkyboxColor = texture(samplerCubeMap, reflect(-V, N)).rgb;
    
    // Reflectance equation
    vec3 Lo = vec3(0.0);
    for (int i = 0; i < ubd.maxLight; i++)
    {
        float shadow = 1.0;
        vec3 L = normalize(ubl[i].position - WorldPos);
        vec3 H = normalize(V + L);
        float distance = length(ubl[i].position - WorldPos);
        float attenuation = 1.0 / (distance * distance);
        vec3 radiance = ubl[i].color * attenuation * ubl[i].range;
        if(ubl[i].status == 0)
        {
           L = -ubl[i].direction;
           H = normalize(V + L);
        }

        // Cook-Torrance BRDF
        float NDF = DistributionGGX(N, H, roughness);        
        float G = GeometrySmith(N, V, L, roughness);      
        vec3 F = fresnelSchlick(max(dot(H, V), 0.0), F0);       

        vec3 kS = F;
        vec3 kD = vec3(1.0) - kS;
        kD *= 1.0 - metallic;      

        vec3 numerator = NDF * G * F;
        float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.0001;
        vec3 specular = (numerator / denominator);  
            
        // Add to outgoing radiance Lo
        float NdotL = max(dot(N, L), 0.0);

        if (ubl[i].status == 0) // DirLight
        {
            if(ubl[i].shadowID >= 0 && ubm[imaterial].castShadow == 1)
            {
                shadow = isInShadow((ubs[ubl[i].shadowID].proj * ubs[ubl[i].shadowID].view) * vec4(WorldPos,1.0),ubl[i].shadowID);
            }
            Lo += (kD * color + specular) * ubl[i].color * ubl[i].range/10.0f * NdotL * shadow;
        }
        else if (ubl[i].status == 1) // PointLight
        {
            Lo += (kD * color / PI + specular) * radiance * NdotL;
        }
        else if (ubl[i].status == 2) // SpotLight
        {
            vec3 lightDir = normalize(ubl[i].direction);
            float spotAngle = radians(ubl[i].spotAngle);
            float spotEffect = dot(lightDir, -L);            

            if(ubl[i].shadowID >= 0 && ubm[imaterial].castShadow == 1)
            {
                shadow = isInShadow((ubs[ubl[i].shadowID].proj * ubs[ubl[i].shadowID].view) * vec4(WorldPos,1.0),ubl[i].shadowID);
            }

            if (spotEffect > cos(spotAngle / 2.0))
            {
                float transitionAngle = radians(4.0);
                float edge0 = cos(spotAngle / 2.0 - transitionAngle);
                float edge1 = cos(spotAngle / 2.0);
                float smoothFactor = smoothstep(edge1, edge0, spotEffect);

                Lo += (kD * color / PI + specular) * radiance * NdotL * pow(smoothFactor, 2.0) * shadow;
            }
        }        
    }    

    color = ambient + Lo*n;
        
    color = mix((color*reflectedSkyboxColor),color,roughness);

    color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0/2.2));  

    outColor = vec4(color, cr.a);
}