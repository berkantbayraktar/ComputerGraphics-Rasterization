#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "hw2_types.h"
#include "hw2_math_ops.h"
#include "hw2_file_ops.h"
#include <iostream>


Camera cameras[100];
int numberOfCameras = 0;

Model models[1000];
int numberOfModels = 0;

Color colors[100000];
int numberOfColors = 0;

Translation translations[1000];
int numberOfTranslations = 0;

Rotation rotations[1000];
int numberOfRotations = 0;

Scaling scalings[1000];
int numberOfScalings = 0;

Vec3 vertices[100000];
int numberOfVertices = 0;

Color backgroundColor;

// backface culling setting, default disabled
int backfaceCullingSetting = 0;

Color **image;


// Helper functions
void translate_triangle(Triangle* triangle,Translation translation);
void rotate_triangle(Triangle * triangle,Rotation rotation);
void scale_triangle(Triangle* triangle,Scaling scaling);

/*
	Initializes image with background color
*/
void initializeImage(Camera cam) {
    int i, j;

    for (i = 0; i < cam.sizeX; i++)
        for (j = 0; j < cam.sizeY; j++) {
            image[i][j].r = backgroundColor.r;
            image[i][j].g = backgroundColor.g;
            image[i][j].b = backgroundColor.b;

        }
}

/*
	Transformations, culling, rasterization are done here.
	You can define helper functions inside this file (rasterizer.cpp) only.
	Using types in "hw2_types.h" and functions in "hw2_math_ops.cpp" will speed you up while working.
*/
void forwardRenderingPipeline(Camera cam) {
    // Traverse all models
    for(int i = 0; i < numberOfModels;i++)
    {   
        // pointer for model, we want to change this model in functions
        // Don't forget te deallocate
        Model* model = new Model();
        model = &models[i];
        // Traverse model's transformations
        for(int j = 0 ; j < model -> numberOfTransformations ;j++)
        {   
            //Traverse model's triangles
            for(int z = 0 ;z < model -> numberOfTriangles;z++)
            {
                // pointer for triangle, we want to change this triangle in functions
                // Don't forget te deallocate
                Triangle* triangle = new Triangle();
                triangle = &(model -> triangles[z]);
                
                // transform id of
                int transform_id = model -> transformationIDs[j];
                //transformation type (rotation,scale,translate)
                char type = model -> transformationTypes[j];
                // If translation
                if(type == 't')
                {   
                    // translation amounts of x,y,z
                    Translation translation = translations[transform_id - 1];
                    translate_triangle(triangle,translation);
                }
                // If rotation
                else if(type == 'r')
                {   
                    //rotation angle, and axis
                    Rotation rotation = rotations[transform_id - 1];
                    rotate_triangle(triangle,rotation);
                }
                // If scaling
                else if(type == 's')
                {   
                    //scale factors
                    Scaling scaling =  scalings[transform_id - 1];
                    scale_triangle(triangle,scaling);
                }
            
            }
        }
    }

}

void translate_triangle(Triangle* triangle,Translation translation)
{

}


void rotate_triangle(Triangle * triangle,Rotation rotation)
{

}



void scale_triangle(Triangle* triangle,Scaling scaling)
{

}



int main(int argc, char **argv) {
    int i, j;

    if (argc < 2) {
        std::cout << "Usage: ./rasterizer <scene file> <camera file>" << std::endl;
        return 1;
    }

    // read camera and scene files
    readSceneFile(argv[1]);
    readCameraFile(argv[2]);

    image = 0;

    for (i = 0; i < numberOfCameras; i++) {

        // allocate memory for image
        if (image) {
			for (j = 0; j < cameras[i].sizeX; j++) {
		        delete image[j];
		    }

			delete[] image;
		}

        image = new Color*[cameras[i].sizeX];

        if (image == NULL) {
            std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
            exit(1);
        }

        for (j = 0; j < cameras[i].sizeX; j++) {
            image[j] = new Color[cameras[i].sizeY];
            if (image[j] == NULL) {
                std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
                exit(1);
            }
        }


        // initialize image with basic values
        initializeImage(cameras[i]);

        // do forward rendering pipeline operations
        forwardRenderingPipeline(cameras[i]);

        // generate PPM file
        writeImageToPPMFile(cameras[i]);

        // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
        // Notice that os_type is not given as 1 (Ubuntu) or 2 (Windows), below call doesn't do conversion.
        // Change os_type to 1 or 2, after being sure that you have ImageMagick installed.
        convertPPMToPNG(cameras[i].outputFileName, 99);
    }

    return 0;

}
