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

    // Convert texCoord to pixel index. For TextureNearest you only need to use the integer part to find it's index. (indices range from 0 to width/height - 1)

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
    // Pixel Location - Same method as NN to find index of top left corner of pixel 

    int pixelx_idx = std::floor(texCoord.x * (image.width));
    int pixely_idx = std::floor(texCoord.y * (image.height));

    pixelx_idx = glm::clamp(pixelx_idx, 0, image.width - 1);
    pixely_idx = glm::clamp(pixely_idx, 0, image.height - 1);

    // Saving offset for interpolation (this will be the alpha and beta values in bilinear interpolation equation). 
    // The offset is calculated by subtracting the integer value from the full decimal value of the coordinate
    // IE: in 4x4 image: (1,1) maps to (4,4). Thus it will have offset(1,1). 

    int offsetX = (texCoord.x * (image.width)) - pixelx_idx; //alpha
    int offsetY = (texCoord.y * (image.height)) - pixely_idx; //beta

    // Indexes of other neighbours. 
    // Make sure to clamp to ensure there is no unexpected behaviour at the borders.

    int x1 = glm::clamp(pixelx_idx + 1, 0, image.width - 1);
    int y1 = glm::clamp(pixely_idx + 1, 0, image.height - 1);

    // Retrieve values of neighbouring texels using indexes calculated using the formula from the lecture

    glm::vec3 texelTL = image.pixels[pixely_idx * image.width + pixelx_idx];
    glm::vec3 texelTR = image.pixels[pixely_idx * image.width + x1];
    glm::vec3 texelBL = image.pixels[y1 * image.width + pixelx_idx];
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

    return interpolatedTexture;
}