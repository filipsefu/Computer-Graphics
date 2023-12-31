#include "interpolate.h"
#include <glm/geometric.hpp>

float findArea(const glm::vec3& v0, const glm::vec3& v1)
{
    //Use equation: .5 * the magnitude of the cross of 2 / 3 vectors which make up the triangle
    return 0.5f * glm::length(glm::cross(v0, v1));
}

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    //Computer Graphics 2023/2024 Lecture 5, slide 81 (barycentric coordinates)
    //Initialize
    float alpha;
    float beta;
    float gamma;

    //Compute vectors
    glm::vec3 v0P = p - v0;
    glm::vec3 v0v1 = v1 - v0;
    glm::vec3 v0v2 = v2 - v0;

    glm::vec3 v1P = p - v1;
    glm::vec3 v1v2 = v2 - v1;

    glm::vec3 v2P = p - v2;
    glm::vec3 v2v0 = v0 - v2;

    //Compute areas
    float totalArea = findArea(v0v1, v0v2);

    float areaAlpha = findArea(v1P, v1v2);
    float areaBeta = findArea(v2P, v2v0);
    float areaGamma = findArea(v0P, v0v1);

    //Prevent division by 0
    if (totalArea == 0) {
        return glm::vec3(0);
    }
    //Compute Barycentric Coordinates
    alpha = areaAlpha / totalArea;
    beta = areaBeta / totalArea;
    gamma = areaGamma / totalArea;

    //Return
    return glm::vec3(alpha, beta, gamma);
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    glm::normalize(n0);
    glm::normalize(n1);
    glm::normalize(n2);
    //bc contains the weights of each normal so simply multiply them accordingly, then add.
    glm::vec3 interpolatedNormal = n0 * bc.x + n1 * bc.y + n2 * bc.z;
    return glm::normalize(interpolatedNormal);
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
    //check interpolateNormal for explanation
    glm::vec2 interpolatedCoord = t0 * bc.x + t1 * bc.y + t2 * bc.z; 
    return interpolatedCoord;
}


