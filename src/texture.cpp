#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image

    // Convert texCoord to pixel index (multiply by image resolution)
    int pixelx_idx = std::floor(texCoord.x * (image.width));
    int pixely_idx = std::floor(texCoord.y * (image.height));

    // Clamp to image indices (texture coordinate (1,1) maps to (4,4) but should map to (3,3)), this will set the maximum allowed value to (3,3).
    pixelx_idx = glm::clamp(pixelx_idx, 0, image.width - 1);
    pixely_idx = glm::clamp(pixely_idx, 0, image.height - 1);

    // Convert (i,j) to index using lecture method
    int pixelIndex = pixely_idx * image.width + pixelx_idx;

    return image.pixels[pixelIndex];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // Pixel Location
    int pixelx_idx = std::floor(texCoord.x * (image.width));
    int pixely_idx = std::floor(texCoord.y * (image.height));

    // Saving offset for interpolation
    int offsetX = (texCoord.x * (image.width)) - std::floor(texCoord.x * (image.width));
    int offsetY = (texCoord.y * (image.height)) - std::floor(texCoord.y * (image.height));

    // Top - left pixel idx
    int tlx_idx = std::floor(texCoord.x * (image.width - 1));
    int tly_idx = std::floor(texCoord.y * (image.height - 1));

    // Other idx (offset for other neighbours)
    // Use clamp to ensure we stay within the image border
    int x1 = glm::clamp(tlx_idx + 1, 0, image.width - 1);
    int y1 = glm::clamp(tly_idx + 1, 0, image.height - 1);

    // Retrieve values of neighbouring texels using formula from lecture

    glm::vec3 texelTL = image.pixels[tly_idx * image.width + tlx_idx];
    glm::vec3 texelTR = image.pixels[tly_idx * image.width + x1];
    glm::vec3 texelBL = image.pixels[y1 * image.width + tlx_idx];
    glm::vec3 texelBR = image.pixels[y1 * image.width + x1];

    // Bilinear Interpolation
    //  ((1-B) * (((1-A) * TopLeft) + ((A) * TopRight))) + (B * ((1-A) * BottomLeft + (A) * BottomRight))

    glm::vec3 interpolatedTexture = ((1.0f - offsetY) * //(1-B)
                                        ((1.0f - offsetX) * texelTL + //(1-A) * TopLeft
                                            ((0.0f + offsetX) * texelTR))) //(A) * TopRight
        + ((0.0f + offsetY) * //(B)
            ((1.0f - offsetX) * texelBL + //(1-A) * BottomLeft
                ((0.0f + offsetX) * texelBR))) //(A) * BottomRight
        ;

    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image

    return interpolatedTexture;
}