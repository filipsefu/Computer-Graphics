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
#include <draw.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
//
// reference:
// https://en.wikipedia.org/wiki/Depth_of_field
// "Fundamentals of Computer Graphics - Ch 13.4.3 Depth of Field"
// //
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int j = 0; j < screen.resolution().y; j++) {
        for (int i = 0; i < screen.resolution().x; i++) {
            uint16_t sampler = screen.indexAt(i, j);
            RenderState state = { .scene = scene, .features = features, .bvh = bvh, .sampler = sampler };
            std::vector<Ray> rays = generatePixelRays(state, camera, { i, j }, screen.resolution());
            std::vector<Ray> rays2;

            glm::vec3 point;
            glm::vec2 noise;

            // for visual debug
            float focalLength = features.focalLength;
            bool visualize = false;
            if (visualize)
                focalLength = glm::sqrt(11.0f);

            for (Ray& r : rays) {
                // focal point
                point = r.origin + r.direction * focalLength;

                for (int k = 0; k < 10; k++) {
                    noise = state.sampler.next_2d() - 0.5f;
                    r.origin += glm::vec3(noise.x * features.aperture, noise.y * features.aperture, 0.0f);
                    r.direction = glm::normalize(-r.origin + point);
                    //res.t = r.t;

                    rays2.push_back(r);
                }
            }
            screen.setPixel(i, j, renderRays(state, rays2));
        }
    }
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!

//Sources that aided my general understanding of motion blur:
//https://en.wikipedia.org/wiki/Motion_blur - specifically “Photography”
//Fundamentals of Computer Graphics fourth edition 13.4.5 (page 349-350)
//Computer Graphics 2023/2024 Lecture 8 - Ray Tracing

//I realized thursday (Nov 2) that my implementation did not take the right direction. However, the correct implementation would require me to change sensitive files like CMAKE and 
//other files that are used by the entire project so I did not want to risk putting my group's project in danger of breaking at the last second. Thus, I am leaving my original implementation.
//My other implementation would consist of giving every ray a sampled time (in generate pixelRays) and altering the intersect functions to move the object based on the bezier curve, then using renderRays to get the average color at each pixel.
//Even though my implementation is not fully correct I will add some visual debug.

// CubicBezierCurve struct holds the positions that define the motion path
// index 0 is start, index 3 is end.

struct CubicBezierCurve {
    std::vector<glm::vec3> keyframes;
};

void interpolateMesh(Mesh& mesh, CubicBezierCurve curve, float t)
{

    //To figure out the formula I used the following video: https://www.youtube.com/watch?v=pnYccz1Ha34&pp=ygUMYmV6aWVyIGN1cnZl.
;
    float tt = t * t;
    float ttt = tt * t;
    float u = 1.0f - t;
    float uu = u * u;
    float uuu = u * u * u;

    // Cubic Bezier interpolation.
    // (1-t)^3 * start + (3((1-t)^2)*t) * kf1 + (3(1-t) * (t**2)) * kf2 + t**3 * end

    for (int i = 0; i < mesh.vertices.size(); i++) {
        glm::vec3 originalPos = mesh.vertices[i].position;
        glm::vec3 interpolatedPos = (uuu * curve.keyframes[0]) + (3 * uu * t * curve.keyframes[1]) + (3 * u * tt * curve.keyframes[2]) + (ttt * curve.keyframes[3]);

        // update vertex, add the translation to the original position
        mesh.vertices[i].position = originalPos + interpolatedPos;
    }
}

void renderFrameWithMotionBlur(Scene& frameScene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen, float time, std::vector<glm::vec3>& colors)
{
    int width = screen.resolution().x;
    int height = screen.resolution().y;

    for (Mesh& mesh : frameScene.meshes) {

        // Create Keyframes for arbitrary scene
        // We need 4 keyframes for cubic bezier curve interpolation
        // We set the first keyframe to be the starting point, thus at 0.0f
        // We set the 4th keyframe to be the endpoint, thus it set to the exposureTime.

        CubicBezierCurve curve = { { glm::vec3(0.0f), glm::vec3(0.1, 0, 0), glm::vec3(0, 0.1, 0), glm::vec3(0, 0, .5) } };

        interpolateMesh(mesh, curve, time);


    }
    
    //Create a BVH for this frame because otherwise it will return color values of the original picture.

    BVH frameBVH = BVH(frameScene, features);

    // Render the frameScene, but store colors in buffer instead of drawing pixels.
    // I copied this section from renderImage, making small necessary changes.

    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
            // Note; we seed the sampler for consistenct behavior across frames
            RenderState state = {
                .scene = frameScene,
                .features = features,
                .bvh = frameBVH,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };
            
            auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
            auto L = renderRays(state, rays);

            //Add color to the buffer

            colors[x + y * width] += L;
        }
    }
}

void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }
    // Set exposureTime and numFrames (will add slider if i have time)
    // exposureTime is the amount of time for the motion
    float exposureTime = 1.0f;

    // numFrames is basically the amount of samples throughout the motion we will take (AKA screencapture of the object at time t)
    int numFrames = 50;

    //Initialize color buffer to store colors.
    std::vector<glm::vec3> colors(screen.resolution().x * screen.resolution().y, glm::vec3(0.0f));

    // Iterate over frames
    for (int frame = 0; frame < numFrames; frame++) {

        // Set new scene for each frame so they can be modified and the colors of their pixels can be computed individually.
        // In other words you create an instance of the scene at certain points in the motion.

        Scene motionBlurredScene = scene;

        // Time at which current frame should be rendered
        // Evenly distributed time intervals. IE: for 5 frames and exposure Time .5, we have: (0,1/8,2/8,3/8,4/8)
        // In other words, divides the exposureTime  interval evenly into (#of frames) values

        float time = frame / static_cast<float>(numFrames - 1) * exposureTime;

        renderFrameWithMotionBlur(motionBlurredScene, bvh, features, camera, screen, time, colors);
    }


    // Now colors contains the total color contributions of all the frames, thus now we extract
    // these color values and set the screen pixels.

    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            glm::vec3 color = colors[x + y * screen.resolution().x];
            glm::vec3 averageColor = color / static_cast<float>(numFrames);
            screen.setPixel(x, y, averageColor);
        }
    }
}
    // TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
//
// reference: "https://brightspace.tudelft.nl/d2l/le/content/595314/viewContent/3520316/View"
// //
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    Screen high(image.resolution());
    Screen horizontal(image.resolution());
    Screen box(image.resolution());
    Screen result(image.resolution());

    float epsilon2 = features.epsilon;

    // take big values
    for (int i = 0; i < image.resolution().x; i++)
    {
        for (int j = 0; j < image.resolution().y; j++) {
            int index = image.indexAt(i, j);
            glm::vec3 color = image.pixels().at(index);

            if (!(color.x > epsilon2 && color.y > epsilon2 && color.z > epsilon2))
                color = glm::vec3(0);

            high.setPixel(i, j, color);
        }
    }

    //filter horizontally
    for (int i = 0; i < image.resolution().x; i++) {
        for (int j = 0; j < image.resolution().y; j++) {
            int index0 = high.indexAt(i, 0);
            int index1 = high.indexAt(i, 0);
            int index2 = high.indexAt(i, j);
            int index3 = high.indexAt(i, image.resolution().y - 1);
            int index4 = high.indexAt(i, image.resolution().y - 1);

            if (j - 2 >= 0)
                index0 = high.indexAt(i, j - 2);
            if (j - 1 >= 0)
                index1 = high.indexAt(i, j - 1);
            if (j + 1 < image.resolution().y)
                index3 = high.indexAt(i, j + 1);
            if (j + 2 < image.resolution().y)
                index4 = high.indexAt(i, j + 2);

            glm::vec3 color = high.pixels().at(index0) + high.pixels().at(index1) * 4.0f + high.pixels().at(index2) * 6.0f +
                high.pixels().at(index3) * 4.0f + high.pixels().at(index4);
            color /= 16;
            horizontal.setPixel(i, j, color);
        }
    }

    // filter vertically
    for (int i = 0; i < image.resolution().x; i++) {
        for (int j = 0; j < image.resolution().y; j++) {
            int index0 = horizontal.indexAt(0, j);
            int index1 = horizontal.indexAt(0, j);
            int index2 = horizontal.indexAt(i, j);
            int index3 = horizontal.indexAt(image.resolution().x - 1, j);
            int index4 = horizontal.indexAt(image.resolution().x - 1, j);

            if (i - 2 >= 0)
                index0 = horizontal.indexAt(i - 2, j);
            if (i - 1 >= 0)
                index1 = horizontal.indexAt(i - 1, j);
            if (i + 1 < image.resolution().x)
                index3 = horizontal.indexAt(i + 1, j);
            if (i + 2 < image.resolution().x)
                index4 = horizontal.indexAt(i + 2, j);

            glm::vec3 color = horizontal.pixels().at(index0) + horizontal.pixels().at(index1) * 4.0f + horizontal.pixels().at(index2) * 6.0f + 
                horizontal.pixels().at(index3) * 4.0f + horizontal.pixels().at(index4);
            color /= 16;
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
//Information regarding the implementation got from the book and the additional document
bool glossyDebug = false;
void setGlossyDebug(bool value)
{
    glossyDebug = value;
}

void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{

    glm::vec3 reflectedRay;
    glm::vec3 pointOfIntersection;
    glm::vec3 glossyAccumulator;
    glm::vec3 orthogonalBasis;
    glm::vec3 orthogonalVector;
    glm::vec3 reflectedRayPrime;
    glm::vec3 renderRayResult;
    glm::vec2 sampler;

    float miscelaneous = 0.0006f;

    reflectedRay = glm::reflect(ray.direction, hitInfo.normal);
    pointOfIntersection = ray.origin + ray.t * ray.direction;

   switch ((reflectedRay.x == 0) * 1 + (reflectedRay.y == 0) * 2) {
    case 1:
        orthogonalVector.x = 1;
        orthogonalVector.y = -reflectedRay.z;
        orthogonalVector.z = reflectedRay.y;
        break;
    case 2:
        orthogonalVector.x = -reflectedRay.z;
        orthogonalVector.y = 1;
        orthogonalVector.z = reflectedRay.x;
        break;
    default:
        orthogonalVector.x = reflectedRay.y;
        orthogonalVector.y = -reflectedRay.x;
        orthogonalVector.z = 0;
    }

    orthogonalBasis = glm::cross(orthogonalVector, reflectedRay);

    orthogonalBasis = glm::normalize(orthogonalBasis);

    for (int i = 0; i < state.features.extra.numGlossySamples; i++) {

        // Two random between
        sampler = state.sampler.next_2d();
        float circlez = glm::sqrt(sampler.y * sampler.y + sampler.x * sampler.x);
        float circleRadius = glm::sqrt(circlez) * hitInfo.material.shininess / 64.0f;
        float formulaAngle = 2 * glm::pi<float>() * sampler.x;

        // Calculation of U and V for the formula
        float u = circleRadius * glm::cos(formulaAngle);
        float v = circleRadius * glm::sin(formulaAngle);

        // Calculation fo the formula
        reflectedRayPrime = reflectedRay + u * orthogonalVector + v * orthogonalBasis;
        reflectedRayPrime = glm::normalize(reflectedRayPrime);

        // If the light is not from behind
        float condition = glm::dot(hitInfo.normal, reflectedRayPrime);
        if (condition > 0) {
            Ray newRay = Ray(pointOfIntersection + miscelaneous * reflectedRayPrime, reflectedRayPrime);
            renderRayResult = renderRay(state, newRay ,rayDepth + 1);
            if (glossyDebug)
                drawRay(newRay, renderRayResult);
            glossyAccumulator = glossyAccumulator + renderRayResult;
        }
    }

    hitColor = hitColor + hitInfo.material.ks * glossyAccumulator / (float)state.features.extra.numGlossySamples;

    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
    // Generate an initial specular ray, and base secondary glossies on this ray
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!

// Information regarding the implementation of this function got from https://en.wikipedia.org/wiki/Cube_mapping
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

        int b = nBuckets* ((centroid - centroidBounds.lower[axis]) / (centroidBounds.upper[axis] - centroidBounds.lower[axis]));

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
                if (count0 == 0)
                b0= buckets[j].bounds;
                count0 += buckets[j].count;
            }

            for (int j = i + 1; j < nBuckets; ++j)
            {
                b1 = Union(b1, buckets[j].bounds);
                if (count1 == 0)
                b1 = buckets[j].bounds;
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
        if (i <=  minCostSplitBucket )
        sum = sum + buckets[i].count;
    }
    if (sum == primitives.size() || sum == 0)
        sum =  splitPrimitivesByMedian(aabb, axis, primitives);
    return sum;
}