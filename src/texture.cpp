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

    //Fundamentals of Computer Graphics fourth edition 11.1 (texture mapping)
     
    //Convert texCoord to index. Multiply by width/height - 1 to account for indices starting at 0.
    //Since images stored upside down, invert y coordinate.
    float xPos = (texCoord.x * (image.width - 1));
    float yPos = (1.0f - texCoord.y) * (image.height - 1);

    // Round to nearest integer (corresponds to index)
    int nearestX = std::round(xPos);
    int nearestY = std::round(yPos);

    // Clamp to image indices to ensure good behaviour near borders.

    nearestX = glm::clamp(nearestX, 0, image.width - 1);
    nearestY = glm::clamp(nearestY, 0, image.height - 1);

    // Convert (i,j) to index 
    // Formula for conversion from Computer Graphics 2023/2024 Lecture 2, slide 27
    int pixelIndex = (nearestY * image.width + nearestX);

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

    float xPos = (texCoord.x * (image.width - 1));
    float yPos = (1.0f - texCoord.y) * (image.height - 1);

    // Use floor instead of round such that you always get the top left texel

    int xTL = std::floor(xPos);
    int yTL = std::floor(yPos);

    // Saving offset for interpolation (this will be the alpha and beta values in bilinear interpolation equation). 
    // The offset is calculated by subtracting the integer value from the full decimal value of the coordinate. 
    // Here I do this by subtracting the idx value of the top left texel (i,j) from the true xPos and yPos

    float offsetX = xPos - xTL; //alpha
    float offsetY = yPos - yTL; //beta

    // Indexes of other neighbours. 
    // Make sure to clamp to ensure there is no unexpected behaviour at the borders.

    int x1 = glm::clamp(xTL + 1, 0, image.width - 1);
    int y1 = glm::clamp(yTL + 1, 0, image.height - 1);

    // Retrieve values of neighbouring texels using indexes calculated using the formula from the lecture

    glm::vec3 texelTL = image.pixels[yTL * image.width + xTL];
    glm::vec3 texelTR = image.pixels[yTL * image.width + x1];
    glm::vec3 texelBL = image.pixels[y1 * image.width + xTL];
    glm::vec3 texelBR = image.pixels[y1 * image.width + x1];

    // Bilinear Interpolation
    //  ((1-B) * (((1-A) * TopLeft) + ((A) * TopRight))) + (B * ((1-A) * BottomLeft + (A) * BottomRight))
    // Formula from Computer Graphics 2023/2024 Lecture 4, slide 76


    //Linear interpolation on upper section
    glm::vec3 interpolationUpper = (1.0f - offsetX) * texelTL + (offsetX) * texelTR;
    //Linear interpolation on lower section
    glm::vec3 interpolationLower = (1.0f - offsetX) * texelBL + (offsetX) * texelBR;
    //Interpolate upper/lower
    glm::vec3 interpolatedTexture = (1.0f - offsetY) * interpolationUpper + (offsetY) * interpolationLower;

    return interpolatedTexture;
}