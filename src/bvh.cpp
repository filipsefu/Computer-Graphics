#include "bvh.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "extra.h"
#include "texture.h"
#include <algorithm>
#include <queue>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>


// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    m_numLevels = 0;
    m_numLeaves = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    m_nodes.emplace_back(); // Create root node
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    buildRecursive(scene, features, primitives, RootIndex, 0);
    // Fill in boilerplate data

    //!-- Optimized and done them when building the BVH --!
    //buildNumLevels();
    //buildNumLeaves();
    m_numLeaves--;

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    AxisAlignedBox aabb;

    // Initialize the AABB corners to the first vertex
    aabb.lower = {
        std::min({ primitive.v0.position.x, primitive.v1.position.x, primitive.v2.position.x }),
        std::min({ primitive.v0.position.y, primitive.v1.position.y, primitive.v2.position.y }),
        std::min({ primitive.v0.position.z, primitive.v1.position.z, primitive.v2.position.z })
    };

    aabb.upper = {
        std::max({ primitive.v0.position.x, primitive.v1.position.x, primitive.v2.position.x }),
        std::max({ primitive.v0.position.y, primitive.v1.position.y, primitive.v2.position.y }),
        std::max({ primitive.v0.position.z, primitive.v1.position.z, primitive.v2.position.z })
    };

    return aabb;
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    AxisAlignedBox globalAABB;
    globalAABB.lower.x = std::numeric_limits<float>::max();
    globalAABB.lower.y = std::numeric_limits<float>::max();
    globalAABB.lower.z = std::numeric_limits<float>::max();

    globalAABB.upper.x = std::numeric_limits<float>::lowest();
    globalAABB.upper.y = std::numeric_limits<float>::lowest();
    globalAABB.upper.z = std::numeric_limits<float>::lowest();

    for (auto prim : primitives) {
        AxisAlignedBox primitiveAABB = computePrimitiveAABB(prim);

        globalAABB.lower.x = std::min(globalAABB.lower.x, primitiveAABB.lower.x);
        globalAABB.lower.y = std::min(globalAABB.lower.y, primitiveAABB.lower.y);
        globalAABB.lower.z = std::min(globalAABB.lower.z, primitiveAABB.lower.z);

        globalAABB.upper.x = std::max(globalAABB.upper.x, primitiveAABB.upper.x);
        globalAABB.upper.y = std::max(globalAABB.upper.y, primitiveAABB.upper.y);
        globalAABB.upper.z = std::max(globalAABB.upper.z, primitiveAABB.upper.z);

    }
    return globalAABB;
}


// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    if (aabb.upper.x - aabb.lower.x >= aabb.upper.y - aabb.lower.y && aabb.upper.x - aabb.lower.x >= aabb.upper.z - aabb.lower.z)
        return 0;
    if (aabb.upper.y - aabb.lower.y >= aabb.upper.x - aabb.lower.x && aabb.upper.y - aabb.lower.y > aabb.upper.z - aabb.lower.z)
        return 1;
    return 2;
}
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{

    float x = (primitive.v0.position.x + primitive.v1.position.x + primitive.v2.position.x) / 3.0f;
    float y = (primitive.v0.position.y + primitive.v1.position.y + primitive.v2.position.y) / 3.0f;
    float z = (primitive.v0.position.z + primitive.v1.position.z + primitive.v2.position.z) / 3.0f;

    glm::vec3 result = glm::vec3(x, y, z);

    return result;
}

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element 
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    std::sort(primitives.begin(), primitives.end(), [axis](const auto& a, const auto& b) {
        return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
    });

   
    if (primitives.size() % 2 == 0) {
        return primitives.size() / 2;
    }
    return primitives.size() / 2 + 1;
}

// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default imlpementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure) {
        std::queue<uint32_t> nodeQueue;
        nodeQueue.push(BVH::RootIndex); // Start with the root node

        while (!nodeQueue.empty()) {
            uint32_t current = nodeQueue.front();
            nodeQueue.pop();

            const BVHInterface::Node& currentNode = nodes[current];

            if (currentNode.isLeaf()) {
                // If the current node leaf -> intersect
                for (uint32_t i = 0; i < currentNode.primitiveCount(); ++i) {
                    const auto& prim = primitives[currentNode.primitiveOffset() + i];
                    const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
                    if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                        updateHitInfo(state, prim, ray, hitInfo);
                        is_hit = true;
                       // drawAABB(nodes[current].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 1.0f);
                    }
                }
            } else {
                // If the current node not a leaf -> check ray intersects  children
                const BVHInterface::Node& leftChild = nodes[currentNode.leftChild()];
                const BVHInterface::Node& rightChild = nodes[currentNode.rightChild()];
                Ray newRay = ray;
                if (intersectRayWithShape(leftChild.aabb, ray)) {
                    ray = newRay;
                    nodeQueue.push(currentNode.leftChild());
                   
                }
                if (intersectRayWithShape(rightChild.aabb, ray)) {
                    ray = newRay;
                    nodeQueue.push(currentNode.rightChild());
                   
                }
            }
        }
    } else {
        // Naive implementation; simply iterates over all primitives
        for (const auto& prim : primitives) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        is_hit |= intersectRayWithShape(sphere, ray, hitInfo);

    return is_hit;
}

// TODO: Standard feature
// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
    Node node;
    node.data[0] = Node::LeafBit;
    node.data[0] |= m_primitives.size();
    node.data[1] = primitives.size();
    node.aabb = aabb;

    for (const auto& primitive : primitives) {
        m_primitives.push_back(primitive);
    }

    return node;
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    Node node;
    node.data[0] = leftChildIndex;
    node.data[1] = rightChildIndex;
    node.aabb = aabb;

    return node;
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex, uint32_t level)
{
    // WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
    // because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.

    // Compute the AABB of the current node.
    // AxisAlignedBox aabb = computeSpanAABB(primitives);

    // As a starting point, we provide an implementation which creates a single leaf, and stores
    // all triangles inside it. You should remove or comment this, and work on your own recursive
    // construction algorithm that implements the following steps. Make sure to reuse the methods
    // you have previously implemented to simplify this process.
    //
    // 1. Determine if the node should be a leaf, when the nr. of triangles is less or equal to 4
    //    (hint; use the `LeafSize` constant)
    // 2. If it is a leaf, fill in the leaf's data, and store its range of triangles in `m_primitives`
    // 3. If it is a node:
    //    3a. Split the range of triangles along the longest axis into left and right subspans,
    //        using either median or SAH-Binning based on the `Features` object
    //    3b. Allocate left/right child nodes
    //        (hint: use `nextNodeIdx()`)
    //    3c. Fill in the current node's data; aabb, left/right child indices
    //    3d. Recursively build left/right child nodes over their respective triangles
    //        (hint; use `std::span::subspan()` to split into left/right ranges)

    // Just configure the current node as a giant leaf for now
    // Compute the AABB of the current node.
    AxisAlignedBox aabb = computeSpanAABB(primitives);

    m_numLevels = std::max(m_numLevels, level + 1);

    const size_t LeafSize = 4; // Adjust this value as needed

    // 1. Determine if the node should be a leaf, when the nr. of triangles is less or equal to 4
    if (primitives.size() <= LeafSize) {
        // 2. If it is a leaf, fill in the leaf's data, and store its range of triangles in `m_primitives`
        {
            m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
            m_numLeaves++;
        }
    } else {
        // 3. If it is a node:
        // 3a. Split the range of triangles along the longest axis into left and right subspans
        uint32_t axis = computeAABBLongestAxis(aabb);
        size_t splitIndex;
        if (!features.extra.enableBvhSahBinning)
            splitIndex = splitPrimitivesByMedian(aabb, axis, primitives);
        else {
            splitIndex = splitPrimitivesBySAHBin(aabb, axis, primitives);
        }

        // 3b. Allocate left/right child nodes
        uint32_t leftChildIndex = nextNodeIdx();
        uint32_t rightChildIndex = nextNodeIdx();

        // 3c. Fill in the current node's data; aabb, left/right child indices
        m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, leftChildIndex, rightChildIndex);

        // 3d. Recursively build left/right child nodes over their respective triangles
        auto leftPrimitives = primitives.subspan(0, splitIndex);
        auto rightPrimitives = primitives.subspan(splitIndex);
        buildRecursive(scene, features, leftPrimitives, leftChildIndex, level + 1);
        buildRecursive(scene, features, rightPrimitives, rightChildIndex, level + 1);
    }
}

// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
    m_numLeaves--;
}

// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int level)
{
    debugDrawLevelHelper(level, 0, 0);
}
void BVH::debugDrawLevelSah(int level)
{
    debugDrawLevelHelperSah(level, 0, 0);
}
double scaleBetweenZeroAndOne(double x, double min, double max)
{
    return (x - min) / (max - min);
}
void BVH::debugDrawLevelHelper(int level, int current_level, int position)
{
    if (position < 0 || position >= m_nodes.size()) {
        return;
    }
    if (m_nodes[position].isLeaf() && current_level != level) {
        return;
    }
    if (level == current_level) {
        drawAABB(m_nodes[position].aabb, DrawMode::Wireframe, glm::vec3(scaleBetweenZeroAndOne(current_level, 0, m_numLevels), 0.5, 1 - scaleBetweenZeroAndOne(current_level, 0, m_numLevels)), 1.0f);
        return;
    }

   if (current_level < level) {
        if (!m_nodes[position].isLeaf()) {

        debugDrawLevelHelper(level, current_level + 1, m_nodes[position].leftChild());
        debugDrawLevelHelper(level, current_level + 1, m_nodes[position].rightChild());
        }
    }


}

// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
    int count = 0;
    debugDrawLeafHelper(leafIndex, 0, count);
}
void BVH::debugDrawLeafSah(int leafIndex)
{
    int count = 0;
    debugDrawLeafHelperSah(leafIndex, 0, count);
}

void BVH::debugDrawLeafHelper(int leafIndex, int position, int& count)
{
    if (position < 0 || position >= m_nodes.size()) {
        return;
    }
    if (m_nodes[position].isLeaf()) {
        if (count == leafIndex) {
            AxisAlignedBox leafAABB = m_nodes[position].aabb;
            drawAABB(leafAABB, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 1.0f);
            return;
        }
        count++;
        return;
    }

    if (count < leafIndex) {
        if (!m_nodes[position].isLeaf()) {
            debugDrawLeafHelper(leafIndex, m_nodes[position].leftChild(), count);
            debugDrawLeafHelper(leafIndex, m_nodes[position].rightChild(), count);
        }
    }
}




void BVH::debugDrawLeafHelperSah(int leafIndex, int position, int& count)
{
    if (position < 0 || position >= m_nodes.size()) {
        return;
    }
    if (m_nodes[position].isLeaf()) {
        if (count == leafIndex) {
            AxisAlignedBox leafAABB = m_nodes[position].aabb;
            drawAABB(leafAABB, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 1.0f);
            return;
        }
        count++;
        return;
    }

    if (count < leafIndex) {
        if (!m_nodes[position].isLeaf()) {
            debugDrawLeafHelperSah(leafIndex, m_nodes[position].leftChild(), count);
            debugDrawLeafHelperSah(leafIndex, m_nodes[position].rightChild(), count);
        }
    }
}

void BVH::debugDrawLevelHelperSah(int level, int current_level, int position)
{
    if (position < 0 || position >= m_nodes.size()) {
        return;
    }
    if (m_nodes[position].isLeaf() && current_level != level) {
        return;
    }
    if (!m_nodes[position].isLeaf() && level == current_level) {
        AxisAlignedBox aabb_left { m_nodes[m_nodes[position].leftChild()].aabb.lower, m_nodes[m_nodes[position].leftChild()].aabb.upper };
        AxisAlignedBox aabb_right { m_nodes[m_nodes[position].rightChild()].aabb.lower, m_nodes[m_nodes[position].rightChild()].aabb.upper };
        drawAABB(aabb_left, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 1.0f);
        drawAABB(aabb_right, DrawMode::Wireframe, glm::vec3(0.0f, 1.0f, 0.0f), 1.0f);
        return;
    }

    if (current_level < level) {
        if (!m_nodes[position].isLeaf()) {
                debugDrawLevelHelperSah(level, current_level + 1, m_nodes[position].leftChild());
                debugDrawLevelHelperSah(level, current_level + 1, m_nodes[position].rightChild());
        }
    }
}