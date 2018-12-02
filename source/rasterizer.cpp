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
bool cull_triangle(Triangle* triangle);
void cam_transform(Triangle *triangle , double M_cam[4][4]);
void per_transform(Triangle *triangle , double M_per[4][4]);
void per_divide(Triangle *triangle);
void vp_transform(Triangle *triangle , double M_vp[3][4]);

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


    double u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,e_x,e_y,e_z,n,f,t,l,b,r,n_x,n_y;
    
    // Camera related points and vectors
    u_x = cam.u.x; u_y = cam.u.y ; u_z = cam.u.z; // u
    v_x = cam.v.x; v_y = cam.v.y ; v_z = cam.v.z; // v
    w_x = cam.w.x; w_y = cam.w.y ; w_z = cam.w.z; // w
    e_x = cam.pos.x; e_y = cam.pos.y ; e_z = cam.pos.z; // e
    n = cam.n; f = cam.f; // near, far
    t = cam.t; b = cam.b; r = cam.r;l = cam.l; // top,right,left,bottom
    n_x = cam.sizeX; n_y = cam.sizeY; // height and width of the camera

    // Camera transformation matrix
    double M_cam[4][4] = { { u_x , u_y , u_z , -(u_x*e_x + u_y*e_y + u_z*e_z) },
                    { v_x , v_y , v_z , -(v_x*e_x + v_y*e_y + v_z*e_z) },
                    { w_x , w_y , w_z , -(w_x*e_x + w_y*e_y + w_z*e_z) },
                    { 0 , 0 , 0 , 1 }};

    // Perspective transformation matrix
    double M_per[4][4] = { { 2*n/(r-l) , 0 , (r+l)/(r-l) , 0 },
                    { 0 , 2*n/(t-b) , (t+b)/(t-b) , 0 },
                    { 0 , 0 , -(f+n)/(f-n) , -(2*f*n) /(f-n)},
                    { 0 , 0 , -1 ,0 }};

    // Viewport transformation matrix
    double M_vp[3][4] = { { n_x/2 , 0 , 0 , (n_x-1)/2 },
                    { 0 , n_y/2 , 0 , (n_y-1)/2 },
                    {0 ,0 , 0.5 , 0.5} };

              
    // Traverse all models
    for(int i = 0; i < numberOfModels;i++)
    {   
        // pointer for model, we want to change this model in functions
        // Don't forget te deallocate
        Model* model = new Model();
        model = &models[i];
        
        // Traverse model's transformations
        for(int j = 0 ;j < model -> numberOfTriangles; j++)
        {   
            // pointer for triangle, we want to change this triangle in functions
            // Don't forget te deallocate
            
            Triangle* triangle = new Triangle();
            triangle = &(model -> triangles[j]);

            //Traverse model's triangles
            for(int k = 0 ; k < model -> numberOfTransformations; k++)
            {
                 
                // transform id of
                int transform_id = model -> transformationIDs[k];
                //transformation type (rotation,scale,translate)
                char type = model -> transformationTypes[k];
                
                // 1) VERTEX PROCESSING (Pipeline 1.step)
                // a ) MODELING TRANSFORMATION
                
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

                // b) CAMERA TRANSFORMATION
                    cam_transform(triangle,M_cam);

                // c) PERSPECTIVE TRANSFORMATION
                    per_transform(triangle,M_per);
                
                // d) PERSPECTIVE DIVIDE 
                    per_divide(triangle);

                // 2) CULLING (Pipeline 3.step)
                
                // Do culling if backfaced culling is enabled
                if(backfaceCullingSetting == 1)
                {
                    bool cull = cull_triangle(triangle);
                    if(cull == 1)
                    {
                        //TO DO : Extunguish triangle

                    }
                }

                // 3) VIEWPORT TRANSFORMATION (Pipeline 4.step)
                    vp_transform(triangle,M_vp);

                // 4) TRIANGLE RASTERIZATION (Pipeline 5.step)


                // FINISH HIM

           
            }
        }
    }
}


void vp_transform(Triangle *triangle , double M_vp[3][4])
{

}

void cam_transform(Triangle *triangle , double M_cam[4][4])
{

}

void per_transform(Triangle *triangle , double M_per[4][4])
{

}

void per_divide(Triangle *triangle)
{

}

bool cull_triangle(Triangle* triangle)
{

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
