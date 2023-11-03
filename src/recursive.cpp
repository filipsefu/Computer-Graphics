#include "recursive.h"
#include "draw.h"
#include "bvh_interface.h"
#include "intersect.h"
#include "extra.h"
#include "light.h"

// This function is provided as-is. You do not have to implement it.
// Given a range of rays, render out all rays and average the result
glm::vec3 renderRays(RenderState& state, std::span<const Ray> rays, int rayDepth)
{
    glm::vec3 L { 0.f };
    for (const auto& ray : rays) {
        L += renderRay(state, ray, rayDepth);
    }
    return L / static_cast<float>(rays.size());
}

// This method is provided as-is. You do not have to implement it.
// Given a camera ray (or secondary ray), tests for a scene intersection, and
// dependent on the results, evaluates the following functions which you must
// implement yourself:
// - `computeLightContribution()` and its submethods
// - `renderRaySpecularComponent()`, `renderRayTransparentComponent()`, `renderRayGlossyComponent()`
glm::vec3 renderRay(RenderState& state, Ray ray, int rayDepth)
{
    // Trace the ray into the scene. If nothing was hit, return early
    HitInfo hitInfo;
    if (!state.bvh.intersect(state, ray, hitInfo)) {
        drawRay(ray, glm::vec3(1, 0, 0));
        return sampleEnvironmentMap(state, ray);
    }

    // Return value: the light along the ray
    // Given an intersection, estimate the contribution of scene lights at this intersection
    glm::vec3 Lo = computeLightContribution(state, ray, hitInfo);

    // Draw an example debug ray for the incident ray (feel free to modify this for yourself)
    drawRay(ray, glm::vec3(1.0f));

    // Given that recursive components are enabled, and we have not exceeded maximum depth,
    // estimate the contribution along these components
    if (rayDepth < 6) {
        bool isReflective = glm::any(glm::notEqual(hitInfo.material.ks, glm::vec3(0.0f)));
        bool isTransparent = hitInfo.material.transparency != 1.f;

        // Default, specular reflections
        if (state.features.enableReflections && !state.features.extra.enableGlossyReflection && isReflective) {
            renderRaySpecularComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Alternative, glossy reflections
        if (state.features.enableReflections && state.features.extra.enableGlossyReflection && isReflective) {
            renderRayGlossyComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Transparency passthrough
        if (state.features.enableTransparency && isTransparent) {
            renderRayTransparentComponent(state, ray, hitInfo, Lo, rayDepth);
        }
    }

    return Lo;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a mirrored ray
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a reflected ray
// This method is unit-tested, so do not change the function signature.
Ray generateReflectionRay(Ray ray, HitInfo hitInfo)
{
    //Normalize ray direction
    glm::vec3 L = glm::normalize(ray.direction);
    glm::vec3 N = glm::normalize(hitInfo.normal);
    
    //Find angle between intersection normal and light direction

    float angleValue = glm::dot(N, L);

    //Calculate reflection vector using the formula L - (2 * dot(normal,light direction) * normal)
    //This formula is from Computer Graphics 2023/2024 Lecture 8, slide 105

    glm::vec3 R = L - 2.0f * angleValue * N;
    glm::vec3 normalizedR = glm::normalize(R);

    //Computer Graphics 2023/2024 Lecture 8, slide 6 (ray representation)

    glm::vec3 intersection = ray.origin + ray.t * ray.direction;

    //Create reflected ray, it's origin will be the intersection point.
    //Add an offset in the reflected ray direction as to prevent immediate self intersection. 
    //This ensures that there will be no infinite recursion and the reflections look good.

    Ray reflectedRay = {intersection + .0001f * normalizedR, R};

    return Ray {reflectedRay};
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a passthrough ray for transparency,
// starting at the intersection point and continuing in the same direction.
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a passthrough ray for transparency
// This method is unit-tested, so do not change the function signature.
Ray generatePassthroughRay(Ray ray, HitInfo hitInfo)
{
    //Use intersection point as the new origin

    glm::vec3 passthroughOrigin = ray.origin + ray.t * ray.direction;
    
    glm::vec3 normalizedDirection = glm::normalize(ray.direction);

    return Ray {passthroughOrigin + .0001f * normalizedDirection, ray.direction};
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a mirrored ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and adding the result times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRaySpecularComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    //Check if reflections disabled

    if (!state.features.enableReflections) {
        //Reflections are disabled
        return;
    }

    //Generate reflection ray

    Ray r = generateReflectionRay(ray, hitInfo);
    
    //Get result, multiply by material.ks, add to hitColor
    glm::vec3 reflectionColor = renderRay(state, r, rayDepth + 1);
    hitColor +=  reflectionColor * hitInfo.material.ks;

    //Ensure hitcolor values stay between 0.0f and 1.0f
    hitColor = glm::clamp(hitColor, 0.0f, 1.0f);
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a passthrough transparent ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and correctly alpha blending the result with the current intersection's hit color
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRayTransparentComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    //Check if tranparency disabled

    if (!state.features.enableTransparency) {
        //Transparency is disabled
        return;
    }

    //Generate passthrough ray
    
    Ray r = generatePassthroughRay(ray, hitInfo);
    
    //Get alpha value of intersection for alpha blending

    float alpha = hitInfo.material.transparency;
    
    //Get passthroughcolor
    glm::vec3 passthroughColor = renderRay(state, r, rayDepth + 1); //* hitInfo.material.kd;

    //Do alpha blending using the alpha value of the intersectionobject.

    glm::vec3 alphaBlend = (1.0f - alpha) * passthroughColor  + (alpha) * hitColor;

    hitColor = alphaBlend;
}