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

// Backface culling setting, default disabled
int backfaceCullingSetting = 0;

Color **image;

// Helper function for matrix multiplication
void mat_mul_with_vec4d(double r[3], double m[3][4], double v[4]);

// Helper functions
void translate_triangle(Triangle* triangle,Translation translation);
void rotate_triangle(Triangle * triangle,Rotation rotation);
void scale_triangle(Triangle* triangle,Scaling scaling);
//void cam_transform(Triangle *triangle , double M_cam[4][4]);
//void per_transform(Triangle *triangle, double M_per_cam[4][4]);
void per_cam_transform(Triangle *triangle, double M_per_cam[4][4]);
//void per_divide(Triangle *triangle);
bool cull_triangle(Triangle* triangle, Vec3 cam_pos);
void vp_transform(Triangle *triangle, double M_vp[3][4]);
void midpoint(Triangle *triangle);
void fill_inside(Triangle *triangle);

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

    // Make perspective and camera transformation together
    double M_per_cam[4][4];

    // M_per_cam is handling both perspective and camera transformation
    multiplyMatrixWithMatrix(M_per_cam,M_per,M_cam);

    // Viewport transformation matrix
    double M_vp[3][4] = { { n_x * 0.5 , 0 , 0 , (n_x-1) * 0.5 },
                    { 0 , n_y * 0.5 , 0 , (n_y-1) * 0.5 },
                    {0 ,0 , 0.5 , 0.5} };

              
    // Traverse all models
    for(int i = 0; i < numberOfModels;i++)
    {   
        // pointer for model, we want to change this model in functions
        Model* model = &(models[i]);
        
        // Traverse model's transformations
        for(int j = 0 ;j < model -> numberOfTriangles; j++)
        {   
            // pointer for triangle, we want to change this triangle in functions
           
            Triangle* triangle = &(model -> triangles[j]);

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
            }    
            // b) CAMERA TRANSFORMATION
            //cam_transform(triangle,M_cam);

            // c) PERSPECTIVE TRANSFORMATION
            //per_transform(triangle,M_per);

            // d) PERSPECTIVE DIVIDE 
            //per_divide(triangle);

            // Make camera ,perspective transformation and perspective 
            // dividing together

            per_cam_transform(triangle,M_per_cam);
                
            // 2) CULLING (Pipeline 3.step)
                
            // Do culling if backfaced culling is enabled
            if(backfaceCullingSetting == 1)
            {
                    

                if(cull_triangle(triangle,cam.pos))
                {
                    // Pass triangle
                    // Don't draw
                    continue;
                }
            }
            // 3) VIEWPORT TRANSFORMATION (Pipeline 4.step)
                vp_transform(triangle,M_vp);

            // 4) TRIANGLE RASTERIZATION (Pipeline 5.step)

                // a ) MIDPOINT ALGORITHM
                midpoint(triangle);
                 
                // IF SOLID FRAME : Use triangle barycentric 
                // coordinates to fill triangle's inside
            if(model -> type == 1)
            {                    
                fill_inside(triangle);
            }

            // FINISH HIM                     
        }
    }
}

void translate_triangle(Triangle* triangle,Translation translation)
{   
    
    double M[4][4] = {{ 1 ,0 , 0 , translation.tx },
                    {0 , 1 , 0 , translation.ty },
                    {0 , 0 , 1 , translation.tz },
                    {0 , 0 , 0 , 1}};

    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    // Make vertices homogenous
    double a_h[4] = { a.x , a.y , a.z , 1 };
    double b_h[4] = { b.x , b.y , b.z , 1 };
    double c_h[4] = { c.x , c.y , c.z , 1 };

    // Create array to store result
    double a_res[4],b_res[4],c_res[4];

    // Multiply each vertex with transformation matrix
    // and store result in dummy arrays
    multiplyMatrixWithVec4d(a_res,M,a_h);
    multiplyMatrixWithVec4d(b_res,M,b_h);
    multiplyMatrixWithVec4d(c_res,M,c_h);

    // crop homogenous parts from result
    a.x = a_res[0]; a.y = a_res[1]; a.z = a_res[2];
    b.x = b_res[0]; b.y = b_res[1]; b.z = b_res[2];
    c.x = c_res[0]; c.y = c_res[1]; c.z = c_res[2];

    // Put transformed vertices 
    vertices[triangle -> vertexIds[0] - 1] = a;
    vertices[triangle -> vertexIds[1] - 1] = b;
    vertices[triangle -> vertexIds[2] - 1] = c;
}
// Rotation is little messy. Do it later
void rotate_triangle(Triangle * triangle,Rotation rotation)
{   
    // Use alternative method for rotating
    Vec3 u = { rotation.ux , rotation.uy , rotation.uz};
    Vec3 v,w;

    // Make smallest one zero while converting to v
    // If ux is smallest
    if(abs(rotation.ux) <= abs(rotation.uy) && abs(rotation.ux) <= abs(rotation.uz))
    {   
        // Make smallest one zero
        v.x = 0;
        // Swap and negate one of them
        v.y = -rotation.uz;
        v.z = rotation.uy;
    }
    // If u_y is smallest
    else if(abs(rotation.uy) <= abs(rotation.ux) && abs(rotation.uy) <= abs(rotation.uz))
    {
        // Make smallest one zero
        v.y = 0;
        // Swap and negate one of them
        v.x = -rotation.uz;
        v.z = rotation.ux;
    }
    // If u_z is smallest
    else if(abs(rotation.uz) <= abs(rotation.uy) && abs(rotation.uz) <= abs(rotation.ux))
    {
        // Make smallest one zero
        v.z = 0;
        // Swap and negate one of them
        v.x = -rotation.uy;
        v.y = rotation.ux;
    }

    // Find w through cross product u x v
    w = crossProductVec3(u,v);

    // normalize v and w
    v = normalizeVec3(v);
    w = normalizeVec3(w);

    // Matrix M
    double M[4][4] = { { u.x , u.y , u.z , 0 },
                    { v.x , v.y , v.z , 0 },
                    { w.x , w.y , w.z , 0 },
                    { 0 , 0 , 0 , 1}};
    // Matrix M^-1 
    double M_reverse[4][4] = { { u.x , v.x , w.x , 0},
                            { u.y , v.y , w.y , 0},
                            { u.z , v.z , w.z , 0},
                            { 0 , 0 , 0 , 1}};

    // Rotation matrix around x
    // I assumed we will rotation around x. I think it is okay
    // It is what the slides says
    double R_x[4][4] = {{ 1 , 0 , 0 , 0},
                    { 0 , cos(rotation.angle) , -sin(rotation.angle) , 0},
                    { 0 , sin(rotation.angle) ,cos(rotation.angle) , 0},
                    { 0 , 0 , 0 , 1}};

    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    // Make vertices homogenous
    double a_h[4] = { a.x , a.y , a.z , 1 };
    double b_h[4] = { b.x , b.y , b.z , 1 };
    double c_h[4] = { c.x , c.y , c.z , 1 };

    // Create array to store result
    double a_res[4],b_res[4],c_res[4];

    // Multiply each vertex with transformation matrix
    // and store result in dummy arrays
    
    // Multiply wih  matrix M
    multiplyMatrixWithVec4d(a_res,M,a_h);
    multiplyMatrixWithVec4d(b_res,M,b_h);
    multiplyMatrixWithVec4d(c_res,M,c_h);

    // Perform actual rotation
    multiplyMatrixWithVec4d(a_res,R_x,a_h);
    multiplyMatrixWithVec4d(b_res,R_x,b_h);
    multiplyMatrixWithVec4d(c_res,R_x,c_h);
  
    // Undo M back
    multiplyMatrixWithVec4d(a_res,M_reverse,a_h);
    multiplyMatrixWithVec4d(b_res,M_reverse,b_h);
    multiplyMatrixWithVec4d(c_res,M_reverse,c_h);
    

    // crop homogenous parts from result
    a.x = a_res[0]; a.y = a_res[1]; a.z = a_res[2];
    b.x = b_res[0]; b.y = b_res[1]; b.z = b_res[2];
    c.x = c_res[0]; c.y = c_res[1]; c.z = c_res[2];

    // Put transformed vertices 
    vertices[triangle -> vertexIds[0] - 1] = a;
    vertices[triangle -> vertexIds[1] - 1] = b;
    vertices[triangle -> vertexIds[2] - 1] = c;

}


void scale_triangle(Triangle* triangle,Scaling scaling)
{   

    // Transform matrix for scaling
    double M[4][4] = {{ scaling.sx , 0 , 0 , 0 },
                    { 0 , scaling.sy , 0 , 0 },
                    { 0 , 0 , scaling.sz , 0},
                    { 0 , 0 , 0 , 1}};

    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    // Make vertices homogenous
    double a_h[4] = { a.x , a.y , a.z , 1 };
    double b_h[4] = { b.x , b.y , b.z , 1 };
    double c_h[4] = { c.x , c.y , c.z , 1 };

    // Create array to store result
    double a_res[4],b_res[4],c_res[4];

    // Multiply each vertex with transformation matrix
    // and store result in dummy arrays
    multiplyMatrixWithVec4d(a_res,M,a_h);
    multiplyMatrixWithVec4d(b_res,M,b_h);
    multiplyMatrixWithVec4d(c_res,M,c_h);

    // crop homogenous parts from result
    a.x = a_res[0]; a.y = a_res[1]; a.z = a_res[2];
    b.x = b_res[0]; b.y = b_res[1]; b.z = b_res[2];
    c.x = c_res[0]; c.y = c_res[1]; c.z = c_res[2];

    // Put transformed vertices 
    vertices[triangle -> vertexIds[0] - 1] = a;
    vertices[triangle -> vertexIds[1] - 1] = b;
    vertices[triangle -> vertexIds[2] - 1] = c;
}

// perspective dividing,perspective transformation,camera transformation altogether 
void per_cam_transform(Triangle *triangle , double M_per_cam[4][4])
{
    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    // Make vertices homogenous
    double a_h[4] = { a.x , a.y , a.z , 1 };
    double b_h[4] = { b.x , b.y , b.z , 1 };
    double c_h[4] = { c.x , c.y , c.z , 1 };

    // Create array to store result
    double a_res[4],b_res[4],c_res[4];

    // Multiply each vertex with  transformation matrix
    // and store result in dummy arrays
    multiplyMatrixWithVec4d(a_res,M_per_cam,a_h);
    multiplyMatrixWithVec4d(b_res,M_per_cam,b_h);
    multiplyMatrixWithVec4d(c_res,M_per_cam,c_h);

    // Store homogenous parts  for perspective dividing
    double divide_a = a_res[3];
    double divide_b = b_res[3];
    double divide_c = c_res[3];


    // crop homogenous parts from result
    // and perform perspective dividing
    a.x = a_res[0] / divide_a ; a.y = a_res[1] / divide_a; a.z = a_res[2] / divide_a;
    b.x = b_res[0] / divide_b; b.y = b_res[1] / divide_b; b.z = b_res[2] / divide_b;
    c.x = c_res[0] / divide_c; c.y = c_res[1] / divide_c; c.z = c_res[2] / divide_c;

    // Put transformed vertices 
    vertices[triangle -> vertexIds[0] - 1] = a;
    vertices[triangle -> vertexIds[1] - 1] = b;
    vertices[triangle -> vertexIds[2] - 1] = c;
}

bool cull_triangle(Triangle* triangle, Vec3 cam_pos)
{
    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    Vec3 mid; // find middle point coordinates of the given triangle
    mid.x = (a.x + b.x + c.x) / 3;
    mid.y = (a.y + b.y + c.y) / 3;
    mid.z = (a.z + b.z + c.z) / 3;
    //calculate surface normal of triangle n = (a-mid) X (b- mid)
    Vec3 surface_normal = crossProductVec3(subtractVec3(a,mid),subtractVec3(b,mid)); 

    if(dotProductVec3(surface_normal,subtractVec3(mid,cam_pos)) > 0 )
        return true;
    else
        return false;
}

void vp_transform(Triangle *triangle , double M_vp[3][4])
{
    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    // Make vertices homogenous
    double a_h[4] = { a.x , a.y , a.z , 1 };
    double b_h[4] = { b.x , b.y , b.z , 1 };
    double c_h[4] = { c.x , c.y , c.z , 1 };

    // Create array to store result
    double a_res[3],b_res[3],c_res[3];

    // Multiply each vertex with  transformation matrix
    // and store result in dummy arrays
    multiplyMatrixWithVec4d(a_res,M_vp,a_h);
    multiplyMatrixWithVec4d(b_res,M_vp,b_h);
    multiplyMatrixWithVec4d(c_res,M_vp,c_h);

    // Assign result to vertex
    a.x = a_res[0]; a.y = a_res[1]; a.z = a_res[2];
    b.x = b_res[0]; b.y = b_res[1]; b.z = b_res[2];
    c.x = c_res[0]; c.y = c_res[1]; c.z = c_res[2];

    // Put transformed vertices 
    vertices[triangle -> vertexIds[0] - 1] = a;
    vertices[triangle -> vertexIds[1] - 1] = b;
    vertices[triangle -> vertexIds[2] - 1] = c;

}

void midpoint(Triangle *triangle)
{   /*        a 
              /\  
             /  \ 
           b/____\c
    */       
   int d,y;
    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    Vec3 triangle_vertices [3] = {a,b,c};
   

    for(int j = 0 ; j < 3 ; j++){
        y = triangle_vertices[j].y;
        d = 2 *  (triangle_vertices[j].y - triangle_vertices[(j+1)%3].y) + triangle_vertices[(j+1)%3].x - triangle_vertices[j].x; 
        double c_r = colors[triangle_vertices[j].colorId-1].r;
        double c_g = colors[triangle_vertices[j].colorId-1].g;
        double c_b = colors[triangle_vertices[j].colorId-1].b;
        double dc_r = (colors[triangle_vertices[(j+1)%3].colorId-1].r - colors[triangle_vertices[j].colorId-1].r) / (triangle_vertices[(j+1)%3].x-triangle_vertices[j].x);
        double dc_g = (colors[triangle_vertices[(j+1)%3].colorId-1].r - colors[triangle_vertices[j].colorId-1].r) / (triangle_vertices[(j+1)%3].x-triangle_vertices[j].x);
        double dc_b = (colors[triangle_vertices[(j+1)%3].colorId-1].r - colors[triangle_vertices[j].colorId-1].r) / (triangle_vertices[(j+1)%3].x-triangle_vertices[j].x);

        for(int i = 0 ; i < triangle_vertices[(j+1)%3].x - triangle_vertices[j].x ; i++){

            Color color ;
            color.r = make_between_0_255(c_r);
            color.g = make_between_0_255(c_g);
            color.b = make_between_0_255(c_b);
            image[(int)triangle_vertices[j].x][(int) triangle_vertices[j].y ] = color;
            if(d < 0){
                y = y+ 1;
                d+= triangle_vertices[j].y- triangle_vertices[(j+1)%3].y + triangle_vertices[(j+1)%3].x - triangle_vertices[j].x;
            }
            else{
                d += triangle_vertices[j].y - triangle_vertices[(j+1)%3].y;
            }
            c_r += dc_r;
            c_g += dc_g;
            c_b += dc_b;
        }
    }
    

}

void fill_inside(Triangle *triangle)
{
    Vec3 a = vertices[triangle -> vertexIds[0] - 1];
    Vec3 b = vertices[triangle -> vertexIds[1] - 1];
    Vec3 c = vertices[triangle -> vertexIds[2] - 1];

    double x_0 = a.x; double y_0 = a.y;
    double x_1 = b.x; double y_1 = b.y;
    double x_2 = c.x; double y_2 = c.y;


    int x_min = std::min( std::min(floor(a.x) , floor(b.x)) , floor(c.x));
    int x_max = std::max( std::max(floor(a.x) , floor(b.x)) , floor(c.x));
    int y_min = std::min( std::min(floor(a.y) , floor(b.y)) , floor(c.y));
    int y_max = std::max( std::max(floor(a.y) , floor(b.y)) , floor(c.y));

    double x0_y1 = x_0 * y_1;
    double x0_y2 = x_0 * y_2;
    double x1_y2 = x_1 * y_2;
    double x2_y1 = x_2 * y_1;
    double x2_y0 = x_2 * y_0;
    double x1_y0 = x_1 * y_0;


    double alpha_triangle = x_0 * (y_1 - y_2) + y_0 * (x_2 - x_1) + x1_y2 - x2_y1;
    double beta_triangle  = x_1 * (y_2 - y_0) + y_1 * (x_0 - x_2) + x2_y0 - x0_y2;
    double gamma_triangle = x_2 * (y_0 - y_1) + y_2 * (x_1 - x_0) + x0_y1 - x1_y0;



    Color c_0 = colors[a.colorId - 1];
    Color c_1 = colors[b.colorId - 1];
    Color c_2 = colors[c.colorId - 1];
    // Initial color
    Color color = { 0 , 0 , 0 };

    double alpha,beta,gamma;
    
    // For better performance -> temporal localization
    // It is not finished
    // alpha = (x_min - 1) * (y_1 - y_2) + (y_min - 1) * (x_2 - x_1) + x1_y2 + x2_y1; 
    // beta =  (x_min - 1) * (y_2 - y_0) + (y_min - 1) * (x_0 - x_2) + x2_y0 + x0_y2;
    // gamma = (x_min - 1) * (y_0 - y_1) + (y_min - 1) * (x_1 - x_0) + x0_y1 + x1_y0;

    for(int y = y_min; y < y_max ; y++)
    {
        for(int x = x_min ; x < x_max ; x++)
        {   
            alpha = (x + 0.5) * (y_1 - y_2) + (y + 0.5) * (x_2 - x_1) + x1_y2 - x2_y1;
            if(alpha / alpha_triangle < 0)
                continue;          
            beta = (x + 0.5) * (y_2 - y_0) + (y + 0.5) * (x_0 - x_2) + x2_y0 - x0_y2;
            if(beta / beta_triangle < 0)
                continue;
            gamma = (x + 0.5) * (y_0 - y_1) + (y + 0.5) * (x_1 - x_0) + x0_y1 - x1_y0;
            if(gamma / gamma_triangle < 0)
                continue;

            // Color multiplying performance can be increased
            color.r = make_between_0_255(alpha * c_0.r + beta * c_1.r + gamma * c_2.r);
            color.g = make_between_0_255(alpha * c_0.g + beta * c_1.g + gamma * c_2.g);
            color.b = make_between_0_255(alpha * c_0.b + beta * c_1.b + gamma * c_2.b);

            image[x][y] = color;
        }
    }
  
}

// Multiplication function for Matrix[3][4] with Vec4d
// It is needed for viewport transformation
void mat_mul_with_vec4d(double r[3], double m[3][4], double v[4]) {
    int i, j;
    double total;
    for (i = 0; i < 3; i++) {
        total = 0;
        for (j = 0; j < 4; j++)
            total += m[i][j] * v[j];
        r[i] = total;
    }
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
