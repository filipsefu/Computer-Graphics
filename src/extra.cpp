#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <algorithm>
#include <cmath>
#include <intersect.h>
#include <texture.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    Screen high(image.resolution());
    Screen box(image.resolution());
    Screen result(image.resolution());
    float epsilon = 0.9f;
    float epsilon2 = 0.7f;

    // take big values
    for (int i = 0; i < image.resolution().x; i++)
    {
        for (int j = 0; j < image.resolution().y; j++) {
            int index = image.indexAt(i, j);
            glm::vec3 color = image.pixels().at(index);

            if (color.x < epsilon2 && color.y < epsilon2 && color.z < epsilon2)
                color = glm::vec3(0);

            /*if (color.x > epsilon)
                color.x = 0.0;
            if (color.y < epsilon)
                color.y = 0.0;
            if (color.z < epsilon)
                color.z = 0.0;*/

            high.setPixel(i, j, color);
        }
    }

    //box filter
    for (int i = 1; i < image.resolution().x - 1; i++) {
        for (int j = 1; j < image.resolution().y - 1; j++) {
            int index1 = high.indexAt(i - 1, j - 1);
            int index2 = high.indexAt(i - 1, j);
            int index3 = high.indexAt(i - 1, j + 1);
            int index4 = high.indexAt(i, j - 1);
            int index5 = high.indexAt(i, j);
            int index6 = high.indexAt(i, j + 1);
            int index7 = high.indexAt(i + 1, j - 1);
            int index8 = high.indexAt(i + 1, j);
            int index9 = high.indexAt(i + 1, j + 1);


            glm::vec3 color = high.pixels().at(index1) + high.pixels().at(index2) + high.pixels().at(index3) + 
                high.pixels().at(index4) + high.pixels().at(index5) + high.pixels().at(index6) +
                high.pixels().at(index7) + high.pixels().at(index8) + high.pixels().at(index9);
            color /= 9;
            box.setPixel(i, j, color);
        }
    }

    //add result to original
    for (int i = 0; i < image.resolution().x; i++) {
        for (int j = 0; j < image.resolution().y; j++) {
            int index = image.indexAt(i, j);
            glm::vec3 color = image.pixels().at(index) + box.pixels().at(index);
            color = glm::vec3(glm::min(1.0f, color.x), glm::min(1.0f, color.y), glm::min(1.0f, color.z));
            result.setPixel(i, j, color);
        }
    }
    image = result;
}




// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!

glm::vec3 calculateOrthogonalVector(const glm::vec3& reflectedRay)
{
    glm::vec3 orthogonalVector;
    if (reflectedRay.x == 0) {
        orthogonalVector.x = 1;
        orthogonalVector.y = -reflectedRay.z;
        orthogonalVector.z = reflectedRay.y;
    } else if (reflectedRay.y == 0) {
        orthogonalVector.x = -reflectedRay.z;
        orthogonalVector.y = 1;
        orthogonalVector.z = reflectedRay.x;
    } else {
        orthogonalVector.x = reflectedRay.y;
        orthogonalVector.y = -reflectedRay.x;
        orthogonalVector.z = 0;
    }
    return orthogonalVector;
}

glm::vec3 calculateReflectedRayPrime(const glm::vec3& reflectedRay, const glm::vec3& orthogonalVector, const glm::vec3& orthogonalBasis, const glm::vec2& sampler, const HitInfo& hitInfo)
{
    float circleRadius = glm::sqrt(sampler.y) * hitInfo.material.shininess / 64.0f;
    float formulaAngle = 2 * glm::pi<float>() * sampler.x;
    float u = circleRadius * glm::cos(formulaAngle);
    float v = circleRadius * glm::sin(formulaAngle);
    glm::vec3 reflectedRayPrime = reflectedRay + u * orthogonalVector + v * orthogonalBasis;
    return glm::normalize(reflectedRayPrime);
}

void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    glm::vec3 reflectedRay = glm::reflect(ray.direction, hitInfo.normal);
    glm::vec3 pointOfIntersection = ray.origin + ray.t * ray.direction;
    glm::vec3 glossyAccumulator;
    glm::vec3 orthogonalVector = calculateOrthogonalVector(reflectedRay);
    glm::vec3 orthogonalBasis = glm::cross(orthogonalVector, reflectedRay);
    orthogonalBasis = glm::normalize(orthogonalBasis);
    glm::vec3 reflectedRayPrime;
    glm::vec3 renderRayResult;
    glm::vec2 sampler;
    float miscelaneous = FLT_EPSILON;

    for (int i = 0; i < state.features.extra.numGlossySamples; i++) {
        sampler = state.sampler.next_2d();
        reflectedRayPrime = calculateReflectedRayPrime(reflectedRay, orthogonalVector, orthogonalBasis, sampler, hitInfo);
        float condition = glm::dot(hitInfo.normal, reflectedRayPrime);
        if (condition > 0) {
            renderRayResult = renderRay(state, Ray(pointOfIntersection + miscelaneous * reflectedRayPrime, reflectedRayPrime), rayDepth + 1);
            glossyAccumulator = glossyAccumulator + renderRayResult;
        }
    }

    hitColor = hitColor + hitInfo.material.ks * glossyAccumulator;
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (!state.features.extra.enableEnvironmentMap) {
        return glm::vec3(0.0f);
    }

    AxisAlignedBox myBox;
    myBox.lower = glm::vec3(-1.0f);
    myBox.upper = glm::vec3(1.0f);

    Ray myRay;
    myRay.direction = glm::vec3(0.0f);
    myRay.origin = ray.direction;
    myRay.t = FLT_MAX;
    
    intersectRayWithShape(myBox, myRay);
    glm::vec3 myPoint = myRay.origin + myRay.direction * myRay.t;

    float myX = std::abs(myPoint.x);
    float myY = std::abs(myPoint.y);
    float myZ = std::abs(myPoint.z);

    int myXPos, myYPos, myZPos;

    if (myPoint.x > 0) {
        myXPos = 1;
    } else {
        myXPos = 0;
    }

    if (myPoint.y > 0) {
        myYPos = 1;
    } else {
        myYPos = 0;
    }

    if (myPoint.z > 0) {
        myZPos = 1;
    } else {
        myZPos = 0;
    }

    float myMaxAxis;
    if (myX >= myY && myX 
        >= myZ) {
        myMaxAxis = myX;
    } else if (myY >= myX && myY >= myZ) {
        myMaxAxis = myY;
    } else {
        myMaxAxis = myZ;
    }

    float myUCoord, myVCoord;
    int myCubeMapFace = 0;

    if (myXPos && myX >= myY && myX >= myZ) {
        myMaxAxis = myX;
        myUCoord = -myPoint.z;
        myVCoord = myPoint.y;
        myCubeMapFace = 0;
    } else if (!myXPos && myX >= myY && myX >= myZ) {
        myMaxAxis = myX;
        myUCoord = myPoint.z;
        myVCoord = myPoint.y;
        myCubeMapFace = 1;
    } else if (myYPos && myY >= myX && myY >= myZ) {
        myMaxAxis = myY;
        myUCoord = myPoint.x;
        myVCoord = -myPoint.z;
        myCubeMapFace = 2;
    } else if (!myYPos && myY >= myX && myY >= myZ) {
        myMaxAxis = myY;
        myUCoord = myPoint.x;
        myVCoord = myPoint.z;
        myCubeMapFace = 3;
    } else if (myZPos && myZ >= myX && myZ >= myY) {
        myMaxAxis = myZ;
        myUCoord = myPoint.x;
        myVCoord = myPoint.y;
        myCubeMapFace = 4;
    } else if (!myZPos && myZ >= myX && myZ >= myY) {
        myMaxAxis = myZ;
        myUCoord = -myPoint.x;
        myVCoord = myPoint.y;
        myCubeMapFace = 5;
    }

    float myUValue = (myUCoord / myMaxAxis + 1.0f);
    float myVValue = (myVCoord / myMaxAxis + 1.0f);
    myUValue = myUValue / 2;
    myVValue = myVValue / 2;


    return sampleTextureNearest(*state.scene.environmentMap[myCubeMapFace], glm::vec2(myUValue, myVValue));
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
struct BucketInfo {
    int count = 0;
    std::vector<BVH::Primitive> prim;
    AxisAlignedBox bounds;
};
AxisAlignedBox Union(AxisAlignedBox a, AxisAlignedBox b)
{
    AxisAlignedBox result;
    result.lower.x = std::min(a.lower.x, b.lower.x);
    result.lower.y = std::min(a.lower.y, b.lower.y);
    result.lower.z = std::min(a.lower.z, b.lower.z);
    result.upper.x = std::max(a.upper.x, b.upper.x);
    result.upper.y = std::max(a.upper.y, b.upper.y);
    result.upper.z = std::max(a.upper.z, b.upper.z);

    return result;
}
constexpr int nBuckets = 16;
//For each primitive in the range, we determine the bucket that its centroid lies in and update the bucket�s bounds to include the primitive�s bounds.
void initializeB(std::span<BVH::Primitive> primitiveInfo,const AxisAlignedBox centroidBounds, uint32_t axis, std::vector<BucketInfo> &buckets)
    {
    for (int i = 0; i < primitiveInfo.size(); ++i) {

        float centroid = computePrimitiveCentroid(primitiveInfo[i])[axis];

        int b = nBuckets * ((centroid - centroidBounds.lower[axis]) / (centroidBounds.upper[axis] - centroidBounds.lower[axis]));

        if (b >= nBuckets)
            b = nBuckets - 1;
        if (b <= 0)
            b = 0;

        buckets[b].count++;
        buckets[b].prim.push_back(primitiveInfo[i]);

        if (buckets[b].count == 1)
            buckets[b].bounds = computePrimitiveAABB(primitiveInfo[i]);
        buckets[b].bounds = Union(buckets[b].bounds, computePrimitiveAABB(primitiveInfo[i]));
    }
}


float SurfaceArea(AxisAlignedBox box) 
    {
    // Compute and return the surface area of the box
    return 2 * ((box.upper.x - box.lower.x) * (box.upper.y - box.lower.y) + (box.upper.x - box.lower.x) * (box.upper.z - box.lower.z) + (box.upper.y - box.lower.y) * (box.upper.z - box.lower.z));
    }
    void ComputeSplitCostsWithProbabilities(std::vector<BucketInfo>& buckets, std::vector<float>& cost, const AxisAlignedBox& parentBox,std::span<BVH::Primitive> primitiveInfo)
{

        for (int i = 0; i < nBuckets - 1; ++i) {
            AxisAlignedBox b0, b1;
            b0.upper.x = -FLT_MAX;
            b0.upper.x = -FLT_MAX;
            b0.upper.x = -FLT_MAX;
            b0.lower.x = FLT_MAX;
            b0.lower.x = FLT_MAX;
            b0.lower.x = FLT_MAX;

            b1.upper.x = -FLT_MAX;
            b1.upper.x = -FLT_MAX;
            b1.upper.x = -FLT_MAX;
            b1.lower.x = FLT_MAX;
            b1.lower.x = FLT_MAX;
            b1.lower.x = FLT_MAX;
        
            int count0 = 0, count1 = 0;
            for (int j = 0; j <= i; ++j)
            {
                b0 = Union(b0, buckets[j].bounds);
                count0 += buckets[j].count;
            }

            for (int j = i + 1; j < nBuckets; ++j)
            {
                b1 = Union(b1, buckets[j].bounds);
                count1 += buckets[j].count;
            }

            if (count0 != 0 && count1 != 0)
                cost[i] = count0 * SurfaceArea(b0) + count1 * SurfaceArea(b1);
            else if (count0 == 0)
                cost[i] =  count1 * SurfaceArea(b1);
            else if (count1 == 0)
                cost[i] = count0 * SurfaceArea(b0);
          }
      
}

size_t FindMinCostSplitBucket(std::vector<float> &cost, float& minCost)
{
    minCost = cost[0];
    size_t minCostSplitBucket = 0;
    for (size_t i = 0; i < nBuckets - 1; ++i) {
            if (cost[i] < minCost) {
            minCost = cost[i];
            minCostSplitBucket = i;
        }
    }
    return minCostSplitBucket;
}
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    std::vector<BucketInfo> buckets(nBuckets);
    for (int i = 0; i < nBuckets; i++)
    {
        buckets[i].bounds.upper.x = -FLT_MAX;
        buckets[i].bounds.upper.x = -FLT_MAX;
        buckets[i].bounds.upper.x = -FLT_MAX;
        buckets[i].bounds.lower.x = FLT_MAX;
        buckets[i].bounds.lower.x = FLT_MAX;
        buckets[i].bounds.lower.x = FLT_MAX;
    }
    // Initialize BucketInfo for SAH partition buckets
    initializeB(primitives, aabb, axis, buckets);

    // Compute costs for splitting after each bucket with probabilities
    std::vector<float> cost(nBuckets - 1);
    ComputeSplitCostsWithProbabilities(buckets, cost, aabb, primitives);

    float minCost;
    size_t minCostSplitBucket = FindMinCostSplitBucket(cost, minCost);

    // Return the best split position
    size_t sum = 0;
    size_t poz = 0;

    for (int i = 0; i < nBuckets; i++)
    {
        for (int j = 0; j < buckets[i].prim.size(); j++)
        {
            primitives[poz] = buckets[i].prim[j];
            poz++;
        }
        if (i <  minCostSplitBucket )
        sum = sum + buckets[i].count;
    }
    if (sum == primitives.size() || sum == 0)
        return splitPrimitivesByMedian(aabb, axis, primitives);
    return sum;
}