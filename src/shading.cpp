#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded linear gradient. Feel free to modify this
    static LinearGradient gradient = {
        .components = {
            { 0.1f, glm::vec3(215.f / 256.f, 210.f / 256.f, 203.f / 256.f) },
            { 0.22f, glm::vec3(250.f / 256.f, 250.f / 256.f, 240.f / 256.f) },
            { 0.5f, glm::vec3(145.f / 256.f, 170.f / 256.f, 175.f / 256.f) },
            { 0.78f, glm::vec3(255.f / 256.f, 250.f / 256.f, 205.f / 256.f) },
            { 0.9f, glm::vec3(170.f / 256.f, 170.f / 256.f, 170.f / 256.f) },
        }
    };

    if (state.features.enableShading) {
        switch (state.features.shadingModel) {
            case ShadingModel::Lambertian:
                return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::Phong:
                return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::BlinnPhong:
                return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::LinearGradient:
                return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a Lambertian diffuse shading, returning the reflected light towards the target.
glm::vec3 computeLambertianModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    glm::vec3 N = normalize(hitInfo.normal);
    glm::vec3 L = normalize(lightDirection);

    if (dot(N, L) <= 0)
        return glm::vec3(0.0f, 0.0f, 0.0f);

    return lightColor * sampleMaterialKd(state, hitInfo) * dot(N, L);
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    glm::vec3 N = normalize(hitInfo.normal);
    glm::vec3 L = normalize(lightDirection);

    if (dot(N, L) <= 0)
        return glm::vec3(0.0f, 0.0f, 0.0f);

    // check REFLECTION direction
    glm::vec3 R = normalize(glm::reflect(L, N));
    glm::vec3 V = normalize(cameraDirection);

    if (dot(V, R) <= 0)
        return glm::vec3(0.0f, 0.0f, 0.0f);

    // which way is the lambertian light direction oriented?
    glm::vec3 diffuse = computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);

    return diffuse + hitInfo.material.ks * lightColor * glm::pow(dot(V, R), hitInfo.material.shininess);
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Blinn-Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBlinnPhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    glm::vec3 N = normalize(hitInfo.normal);
    glm::vec3 L = normalize(lightDirection);

    if (dot(N, L) <= 0)
        return glm::vec3(0.0f, 0.0f, 0.0f);

    glm::vec3 V = normalize(cameraDirection);
    glm::vec3 H = normalize(L + V);

    if (dot(N, H) <= 0)
        return glm::vec3(0.0f, 0.0f, 0.0f);

    glm::vec3 diffuse = computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);

    return diffuse + hitInfo.material.ks * lightColor * glm::pow(dot(N, H), hitInfo.material.shininess);
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{    
    if (ti <= -1)
        return this->components.front().color;

    if (ti >= 1)
        return this->components.back().color;

    for (int i = 0; i < this->components.size() - 1; i++) {
        if (ti >= this->components.at(i).t && ti <= this->components.at(i + 1).t)
            //shoudl we resize ti to be in [0,1] ?
                return this->components.at(i).color + (ti - this->components.at(i).t) / 2 * this->components.at(i + 1).color;
    }

    return glm::vec3(0.0f);
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);

    glm::vec3 N = normalize(hitInfo.normal);
    glm::vec3 L = normalize(lightDirection);

    if (dot(N, L) <= 0)
        return glm::vec3(0.0f, 0.0f, 0.0f);

    float ti = hitInfo.barycentricCoord.x;
    if (ti == 0)
        ti = hitInfo.barycentricCoord.y;
    if (ti == 0)
        ti = hitInfo.barycentricCoord.z;

    return gradient.sample(ti) * lightColor * dot(N, L);
}