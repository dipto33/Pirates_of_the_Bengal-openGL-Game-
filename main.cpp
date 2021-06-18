#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <windows.h>
#include <iostream>
#include <sstream>
#include "BmpLoader.h"

using namespace std;
int windowWidth = 800;
int windowHeight = 600;
const int rat = 1 * windowWidth / windowHeight;
//int option = 1;
//bool flagScale = false;
float sz = 0.0f,sharkX=2.0f,bx=0.0f,bz=0.0f,bombY=2.0f,kk=0,gz=-15;
float _angle = 0.0f;
GLfloat cx = 1, cy = 0, cz = 20, ex = 0,  ey = 0,ez = -5;
bool light1 = false;
bool light2 = false;
bool light3 = false;

bool ambientM = true;
bool specularM = true;
bool diffuseM = true;

bool ambient = true;
bool specular = true;
bool diffuse = true;
unsigned int ID;

const double PI = 3.14159265389;
bool sc=false;


/* GLUT callback Handlers */


int anglex= 0, angley = 0, anglez = 0;          //rotation angles
int window;
int wired=0;
int shcpt=1;
int animat = 0;
const int L=13;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 20;


GLfloat ctrlpoints[L+1][3] =
{
    {0.02, 1.5, 0.0},
    { 0.2, 1.5, 0.0},{ 0.5, 1.5, 0.0},
    {1.0, 1.5, 0.0}, {1.4, 1.4, 0.0},
    {1.8, 1.4, 0.0},{2.2, 1.4, 0.0},
    {2.6, 1.5, 0.0},{3.0, 1.4, 0.0},
    {3.4, 1.4, 0.0},{3.8, 1.5, 0.0},
    {4.2, 1.4, 0.0},{4.8, 1.4, 0.0},
    {5.4, 1.5, 0.0}

};

float wcsClkDn[3],wcsClkUp[3];
///////////////////////////////
class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
} clkpt[2];
int flag=0,ct=0;
GLint viewport[4]; //var to hold the viewport info
GLdouble modelview[16]; //var to hold the modelview info
GLdouble projection[16]; //var to hold the projection matrix info

//////////////////////////
void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);
void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);
///////////////////////////

void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
    GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates

    //glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;

    //get the world coordinates from the screen coordinates
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;


}
void processMouse(int button, int state, int x, int y)
{
    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN)
    {
        if(flag!=1)
        {
            flag=1;
            clkpt[0].x=x;
            clkpt[0].y=y;
        }


        scsToWcs(clkpt[0].x,clkpt[0].y,wcsClkDn);
        cout<<"\nD: "<<x<<" "<<y<<" wcs: "<<wcsClkDn[0]<<" "<<wcsClkDn[1];
    }
    else if(button==GLUT_LEFT_BUTTON && state==GLUT_UP)
    {
        if (flag==1)
        {
            clkpt[1].x=x;
            clkpt[1].y=y;
            flag=0;
        }
        float wcs[3];
        scsToWcs(clkpt[1].x,clkpt[1].y,wcsClkUp);
        cout<<"\nU: "<<x<<" "<<y<<" wcs: "<<wcsClkUp[0]<<" "<<wcsClkUp[1];

        clikd=!clikd;
    }
}

//control points
long long nCr(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

static void getNormal3p(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat x3, GLfloat y3, GLfloat z3) {
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2 - x1;
    Uy = y2 - y1;
    Uz = z2 - z1;

    Vx = x3 - x1;
    Vy = y3 - y1;
    Vz = z3 - z1;

    Nx = Uy * Vz - Uz * Vy;
    Ny = Uz * Vx - Ux * Vz;
    Nz = Ux * Vy - Uy * Vx;

    glNormal3f(Nx, Ny, Nz);
}

void bottleBezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                getNormal3p(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}

void showControlPoints()
{
    glPointSize(5.0);
    //glColor3f(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i <=L; i++)
        glVertex3fv(&ctrlpoints[i][0]);
    glEnd();
}
void face(GLfloat A[], GLfloat B[], GLfloat C[], GLfloat D[]) {

    glBegin(GL_POLYGON);
    glVertex3fv(A);
    glTexCoord2f(1,1);
    glVertex3fv(B);
    glTexCoord2f(1,0);
    glVertex3fv(C);
    glTexCoord2f(0,0);
    glVertex3fv(D);
    glTexCoord2f(0,1);
    glEnd();

}
void resize(int windowWidth, int windowHeight) {
    glViewport(100, 20, windowHeight * rat, windowHeight);
}
void doMaterial(GLfloat x, GLfloat y, GLfloat z) {
    GLfloat no_mat[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat mat_ambient[] = {x /3, y, z, 1.0};
    GLfloat mat_diffuse[] = {x, y, z, 1.0};
    GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat mat_shininess[] = {60};
    GLfloat globalAmbient[] ={0.1,0.1,0.1,1.0};
glLightModelfv(GL_LIGHT_MODEL_AMBIENT,globalAmbient);
    if (ambientM) {
        glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    } else {
        glMaterialfv(GL_FRONT, GL_AMBIENT, no_mat);
    }
    if (diffuseM) {
        glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    } else {
        glMaterialfv(GL_FRONT, GL_DIFFUSE, no_mat);
    }
    if (specularM) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    } else {
        glMaterialfv(GL_FRONT, GL_SPECULAR, no_mat);
    }
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
}

void light() {
    GLfloat no_light[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat l_amb[] = {.01, .01, .01, 1};
    GLfloat l_dif[] = {.8, .8, .8, 1};
    GLfloat l_spec[] = {.5, .5, .5, 1};

    GLfloat pos[] = {-3,0, -sz, 1};
    GLfloat pos1[] = {3, 0, -sz, 1};
    GLfloat pos2[] = {0, 2, -sz, 1};
    if (!ambient) {
        l_amb[0] = 0;
        l_amb[1] = 0;
        l_amb[2] = 0;
    }
    if (!diffuse) {
        l_dif[0] = 0;
        l_dif[1] = 0;
        l_dif[2] = 0;
    }
    if (!specular) {
        l_spec[0] = 0;
        l_spec[1] = 0;
        l_spec[2] = 0;
    }
    if (!light1) {
        glDisable(GL_LIGHT0);
    }
    if (light1) {
        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_AMBIENT, l_amb);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, l_dif);
        glLightfv(GL_LIGHT0, GL_SPECULAR, l_spec);
        glLightfv(GL_LIGHT0, GL_POSITION, pos);
    }
    if (!light2) {
        glDisable(GL_LIGHT1);
    }
    if (light2) {
        glEnable(GL_LIGHT1);
        glLightfv(GL_LIGHT1, GL_AMBIENT, l_amb);
        glLightfv(GL_LIGHT1, GL_DIFFUSE, l_dif);
        glLightfv(GL_LIGHT1, GL_SPECULAR, l_spec);
        glLightfv(GL_LIGHT1, GL_POSITION, pos1);
    }
    if (light3) {
        glEnable(GL_LIGHT2);
        glLightfv(GL_LIGHT2, GL_AMBIENT, l_amb);
        glLightfv(GL_LIGHT2, GL_DIFFUSE, l_dif);
        glLightfv(GL_LIGHT2, GL_SPECULAR, l_spec);
        glLightfv(GL_LIGHT2, GL_POSITION, pos2);
        GLfloat spot_direction[] = {0.0, -1.0, 0.0};
        glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, spot_direction);
        glLightf(GL_LIGHT2, GL_SPOT_CUTOFF, 30.0);
    }
    if (!light3) {
        glDisable(GL_LIGHT2);
    }
}
GLfloat V[8][3] =
    {
        {-1, 0.5, .5},
        {1, 0.5, .5},
        {1, -0.5, .5},
        {-1, -0.5, .5},

        {-1, 0.5, -.5},
        {1, 0.5, -.5},
        {1, -0.5, -.5},
        {-1, -0.5, -.5}

};
GLfloat c[8][3] =
    {
        {0, 1.5, 1.5},
        {2, 1.5, 1.5},
        {1, .5, 1.5},
        {0, 0.5, 1.5},

        {0, 1.5, .5},
        {2, 1.5, .5},
        {2, 0.5, .5},
        {0, 0.5, .5}

};
static GLfloat v_pyramid[5][3] =
{
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 0.0, 1.0},
    {1.0, 0.0, 0.0},
    {0.5, 2.0, 0.5}
};

static GLubyte p_Indices[4][3] =
{
    {4, 1, 2},
    {4, 2, 3},
    {4, 3, 0},
    {4, 0, 1}
};
static GLubyte quadIndices[1][4] =
{
    {0, 3, 2, 1}
};

void drawpyramid()
{
    //glColor3f(.5,0.50,0.5);
    glBegin(GL_TRIANGLES);
    for (GLint i = 0; i <4; i++)
    {
        //glColor3f(colors[i][0],colors[i][1],colors[i][2]);
        /* getNormal3p(v_pyramid[p_Indices[i][0]][0], v_pyramid[p_Indices[i][0]][1], v_pyramid[p_Indices[i][0]][2],
                     v_pyramid[p_Indices[i][1]][0], v_pyramid[p_Indices[i][1]][1], v_pyramid[p_Indices[i][1]][2],
                     v_pyramid[p_Indices[i][2]][0], v_pyramid[p_Indices[i][2]][1], v_pyramid[p_Indices[i][2]][2]);*/

        glVertex3fv(&v_pyramid[p_Indices[i][0]][0]);glTexCoord2f(0,0);
        glVertex3fv(&v_pyramid[p_Indices[i][1]][0]);glTexCoord2f(1,0);
        glVertex3fv(&v_pyramid[p_Indices[i][2]][0]);glTexCoord2f(1,1);
    }
    glEnd();

    glBegin(GL_QUADS);
    for (GLint i = 0; i <1; i++)
    {
        //glColor3f(colors[4][0],colors[4][1],colors[4][2]);
        /* getNormal3p(v_pyramid[quadIndices[i][0]][0], v_pyramid[quadIndices[i][0]][1], v_pyramid[quadIndices[i][0]][2],
                     v_pyramid[quadIndices[i][1]][0], v_pyramid[quadIndices[i][1]][1], v_pyramid[quadIndices[i][1]][2],
                     v_pyramid[quadIndices[i][2]][0], v_pyramid[quadIndices[i][2]][1], v_pyramid[quadIndices[i][2]][2]);*/

        glVertex3fv(&v_pyramid[quadIndices[i][0]][0]);glTexCoord2f(1,1);
        glVertex3fv(&v_pyramid[quadIndices[i][1]][0]);glTexCoord2f(1,0);
        glVertex3fv(&v_pyramid[quadIndices[i][2]][0]);glTexCoord2f(0,0);
        glVertex3fv(&v_pyramid[quadIndices[i][3]][0]);glTexCoord2f(0,1);
    }
    glEnd();


}


void cube2(GLfloat V0[], GLfloat V1[], GLfloat V2[], GLfloat V3[], GLfloat V4[], GLfloat V5[], GLfloat V6[], GLfloat V7[]) {
    ///front
    doMaterial(0.545, 0.271, 0.075);
    getNormal3p(V0[0], V0[1], V0[2], V3[0], V3[1], V3[2], V2[0], V2[1], V2[2]);
    face(V0, V1, V2, V3);
    ///back
    //glColor3f(1.0f, 0.5f, 1.0f);
    getNormal3p(V4[0], V4[1], V4[2], V5[0], V5[1], V3[2], V6[0], V6[1], V6[2]);
    face(V4, V5, V6, V7);
    ///left
    doMaterial(0.502, 0.000, 0.000);
    getNormal3p(V0[0], V0[1], V0[2], V4[0], V4[1], V4[2], V7[0], V7[1], V7[2]);
    face(V0, V4, V7, V3);
    ///right
    //glColor3f(1.0f, 1.0f, 0.0f);
    getNormal3p(V1[0], V1[1], V1[2], V2[0], V2[1], V2[2], V6[0], V6[1], V6[2]);
    face(V1, V5, V6, V2);
    ///bottom
    doMaterial(0.502, 0.000, 0.000);
    getNormal3p(V3[0], V3[1], V3[2], V7[0], V7[1], V7[2], V6[0], V6[1], V6[2]);
    face(V3, V2, V6, V7);
    /// top
    //glColor3f(1.0f, 0.0f, 1.0f);
    getNormal3p(V0[0], V0[1], V0[2], V1[0], V1[1], V1[2], V5[0], V5[1], V5[2]);
    face(V0, V1, V5, V4);
}
void cube3(GLfloat V0[], GLfloat V1[], GLfloat V2[], GLfloat V3[], GLfloat V4[], GLfloat V5[], GLfloat V6[], GLfloat V7[], GLfloat x, GLfloat y, GLfloat z) {
    ///top
    doMaterial(x, y, z);
    getNormal3p(V0[0], V0[1], V0[2], V3[0], V3[1], V3[2], V2[0], V2[1], V2[2]);
    face(V0, V1, V2, V3);
    /// front
    //glColor3f(1.0f, 0.5f, 1.0f);
    getNormal3p(V4[0], V4[1], V4[2], V5[0], V5[1], V3[2], V6[0], V6[1], V6[2]);
    face(V4, V5, V6, V7);
    ///right
    //glColor3f(1.0f, 1.0f, 1.0f);
    getNormal3p(V0[0], V0[1], V0[2], V4[0], V4[1], V4[2], V7[0], V7[1], V7[2]);
    face(V0, V4, V7, V3);
    ///back
    //glColor3f(1.0f, 1.0f, 0.0f);
    getNormal3p(V1[0], V1[1], V1[2], V2[0], V2[1], V2[2], V6[0], V6[1], V6[2]);
    face(V1, V5, V6, V2);
    ///left
    //glColor3f(0.0f, 1.0f, 1.0f);
    getNormal3p(V3[0], V3[1], V3[2], V7[0], V7[1], V7[2], V5[0], V5[1], V5[2]);
    face(V3, V2, V6, V7);
    /// bottom
    //glColor3f(1.0f, 1.0f, 1.0f);
    getNormal3p(V0[0], V0[1], V0[2], V1[0], V1[1], V1[2], V5[0], V5[1], V5[2]);
    face(V0, V1, V5, V4);
}
void boat()
{
///Boat_Body
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,1);
    glRotatef(90,0,1,0);
    glPushMatrix();
    glRotatef(90,0,1,0);
    glScalef(.2, .3, 1.5);
    glTranslatef(0,-3,0);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glPopMatrix();
///Boat_Corner
    glPushMatrix();
    glTranslatef(.6,-.8,-.2);
    glRotatef(50,0,0,-1);
    glScalef(.2,.4,.4);
    drawpyramid();
    glPopMatrix();
///Boat_Corner
    glPushMatrix();
    glTranslatef(-.7,-.9,-.2);
    glRotatef(50,0,0,1);
    glScalef(.2,.4,.4);
    drawpyramid();
    glPopMatrix();
///Handle
    glPushMatrix();
    glTranslatef(0,-.6,-.1);
    glScalef(.05,.4,.05);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Pal
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,3);
    glTranslatef(0,0,-.1);
    glRotatef(10,0,-1,0);
    glScalef(.5,1,.05);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7],0.902, 0.902, 0.980);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}
void cloud()
{
    glPushMatrix();
    glPushMatrix();
    glutSolidSphere (.2, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //lColor3f(1,1,0);
    glTranslatef(-.15,.3,0);
    glutSolidSphere (.2, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.1,.5,0);
    glutSolidSphere (.2, 20, 16);
    glPopMatrix();
    //glColor3f(1,1,0);
    glTranslatef(0.3,.4,0);
    glutSolidSphere (.2, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.3,.1,0);
    glutSolidSphere (.2, 20, 16);
    glPopMatrix();
     glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(0,.2,0);
    glutSolidSphere (.2, 20, 16);
    glPopMatrix();
    glPopMatrix();
}
void  tree_body()
{
///Main_Body
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,8);
    glTranslatef(0,-1,.4);
    glScalef(.1,1.5,.1);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Branch
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,8);
    glTranslatef(-.1,-.1,.4);
    glRotatef(30,0,0,1);
    glScalef(.05,.4,.1);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Branch
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,8);
    glTranslatef(.1,-.1,.4);
    glRotatef(30,0,0,-1);
    glScalef(.05,.4,.1);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}
void leaves()
{
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,5);
    //glColor3f(1,1,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //lColor3f(1,1,0);
    glTranslatef(-.15,.2,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(-.2,.35,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(-.15,.5,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(-.1,.65,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(0.1,.7,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.3,.7,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.5,.65,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.6,.45,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.6,.25,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.5,.1,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.2,0,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.3,0,0);
    glutSolidSphere (.15, 20, 16);
    glPopMatrix();
    glPushMatrix();
    //glColor3f(1,1,0);
    glTranslatef(.23,.34,0);
    glutSolidSphere (.32, 20, 16);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}
void  tree()
{
    glPushMatrix();
    glTranslatef(-.2,.05,.35);
    leaves();
    glPopMatrix();
    tree_body();
}

void windmill()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 5000.0;
    const double a = t*90.0;

    glPushMatrix();
    glTranslatef(-2.5,0,-2);
    glRotatef(90,0,1,0);
    glPushMatrix();

    if(animat)
        glRotated(a,0,0,1);
        glScalef(.1,.6,.1);

    glRotatef( anglex, 1.0, 0.0, 0.0);
    glRotatef( angley, 0.0, 1.0, 0.0);         	//rotate about y-axis
    glRotatef( anglez, 0.0, 0.0, 1.0);

    glTranslated(0,-3,0);
    glRotatef( 90, 0.0, 0.0, 1.0);

    glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info

    //matColor(0.9,0.5,0.1,20);   // front face color
    //matColor(0.0,0.5,0.8,20,1);  // back face color


    bottleBezier();
    glPopMatrix();

     glPushMatrix();
     glTranslatef(-.2,1,.2);
   glRotatef(90,1,0,0);
   glRotatef(_angle, 0.0f, -1.0f, 0.0f);
    glPushMatrix();
    glScalef(.5, .05, .15);
    glTranslatef(-1.2, 0, 0);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.412, 0.412, 0.412);
    glPopMatrix();
    //blade front
    glPushMatrix();
    glScalef(.09, .05, 1);
    glTranslatef(0, 0, .5);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.412, 0.412, 0.412);
    glPopMatrix();
    //blade right
    glPushMatrix();
    glScalef(.5, .05, .15);
    glTranslatef(1.2, 0, 0);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.412, 0.412, 0.412);
    glPopMatrix();
    //blade back
    glPushMatrix();
    glScalef(.09, .05, 1);
    glTranslatef(0, 0, -.5);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.412, 0.412, 0.412);
    glPopMatrix();
    //center
    glPushMatrix();
    //glTranslatef(0.01,-.3,0);
    //glRotatef(65,1,0,0);
    glScalef(.1, .1, .15);
    //glRotatef(_angle, 0.0, 1.0, 0);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 1, 1, 1);
    glPopMatrix();
    glPopMatrix();
     glPopMatrix();


}
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(cx, cy, 5-sz, ex, ey, -sz, 0, 1, 0);
    light();


    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //glViewport(0,0,800,800);
///Sea
    for(int i=0;i<50;i++)
    {
    glPushMatrix();
    glTranslatef(0,0,-(i*15));
///Right gift-box
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,6);
    glScalef(.2,.3,.3);
    glTranslatef(8,-3.5,-15);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glPopMatrix();
///middle gift box
    glPushMatrix();
    glScalef(.2,.3,.3);
    glTranslatef(-4,-4,5);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glPopMatrix();
///Left gift box
    glPushMatrix();
    glScalef(.2,.3,.3);
    glTranslatef(-8,-4,25);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Main Sea
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,2);
    glTranslatef(0, -1.5, 0);
    glScalef(2.5, .05, 15);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    }
///Cloud
for(int i=0;i<50;i++)
{
    glPushMatrix();
    glTranslatef(1.5,2,-(i*10));
    cloud();
    glPopMatrix();
    glPushMatrix();
    glTranslatef(-1.5,2,-(i*15));
    cloud();
    glPopMatrix();
         glPushMatrix();
    //glRotatef(90,0,1,0);
glTranslatef(0,0,-(i*15));
    windmill();
    glPopMatrix();

      glPushMatrix();
   glTranslatef(1,-.4,-(30*i));
  boat();
    glPopMatrix();

glPushMatrix();
   glTranslatef(-1,-.4,-(20*i));
  boat();
    glPopMatrix();


///Tree
    glPushMatrix();
    glTranslatef(-2.5,0,-(i*15));
    tree();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-2.5,0,-(i*5));
    tree();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,0,-(i*15));
    tree();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,0,-(i*10));
    tree();
    glPopMatrix();
///Shark
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,4);
    glTranslatef(sharkX,-1.8,-(i*30));
    //cout<<sharkX<<endl;
    glScalef(.3,.4,.3);
    drawpyramid();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

glPushMatrix();
    glTranslatef(0,4,-10);
    doMaterial(1,1,0);
    glTranslatef(0,0,-sz);
    glutSolidSphere (.4, 20, 16);
    glPopMatrix();
///Boat
    glPushMatrix();
    glTranslatef(bx,0,bz);
    glTranslatef(0,-.2,-sz);
    //cout<<bx<<" " <<bz<<endl;;
    boat();
    glPopMatrix();
///Bomb
    for(int i=0;i<50;i++)
    {
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,9);
    glTranslatef(0,bombY,-8);
    //cout<<bombY<<endl;
    glTranslatef(0,0,-(i*10));
    glScalef(.3,.4,.3);
    //glColor3f(1,0,0);
    glutSolidSphere (.3, 20, 16);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,bombY,0);
    glTranslatef(1,0,-(i*15));
    glScalef(.3,.4,.3);
    //glColor3f(1,0,0);
    glutSolidSphere (.3, 20, 16);
    glPopMatrix();
    }
///Left_Sand
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,5);
    glTranslatef(0,0,-sz);
    glTranslatef(-2.5,-1.5,-7);
    glScalef(.2,.3,20);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Right_Sand
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,5);
    glTranslatef(0,0,-sz);
    glTranslatef(2.5,-1.5,-7);
    glScalef(.2,.3,20);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Front Sky
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,7);
    glTranslatef(0,0,-sz);
    glTranslatef(0,8,-15);
    glScalef(5,20,1);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Left Sky
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,7);
    glTranslatef(0,0,-sz);
    glTranslatef(-4,8,-13);
    glScalef(.3,20,30);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
///Right Sky
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,7);
    glTranslatef(0,0,-sz);
    glTranslatef(3.5,8,-13);
    glScalef(.3,20,30);
    cube3(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], 0.902, 0.902, 0.980);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}
void myKeyboardFunc(unsigned char key, int x, int y) {
    switch (key) {
        case 'x':
            cx = cx - 0.5;
            break;
        case 'X':
            cx = cx + 0.5;
            break;
        case 'y':
            cy = cy - 0.5;
            break;
        case 'Y':
            cy = cy + 0.5;
            break;
        case '-':
            cz = cz + .5;
            break;
        case '+':
            cz = cz - .5;
            break;
        case 'L':
            ex = ex - 0.5;
            break;
        case 'R':
            ex = ex + 0.5;
            break;
        case 'D':
            ey = ey - 0.5;
            break;
        case 'U':
            ey = ey + 0.5;
            break;
        case 'b':
            if(bx<=1.5)
            {
                bx = bx + .3;
            }
            break;
        case 'v':
            if(bx>=-1.5)
            {
                bx = bx - .3;
            }
            break;
         case 'n':
            bz = bz - .5;
            break;
        case 'm':
            bz = bz + .5;
            break;
        case '1':
            light1 = !light1;
            break;
        case '2':
            light2 = !light2;
            break;
        case '3':
            light3 = !light3;
            break;
        case 'd':
            diffuse = !diffuse;
            break;
        case 's':
            specular = !specular;
            break;
        case 'a':
            ambient = !ambient;
            break;
        case 'f':
            diffuseM = !diffuseM;
            break;
        case 'g':
            specularM = !specularM;
            break;
        case 'h':
            ambientM = !ambientM;
            break;
    }
    printf("Ambient:%d Specular:%d Diffusion:%d\n" , ambient , specular , diffuse);
    //printf("%d %d %d\n" , ambientM , specularM , diffuseM);
    glutPostRedisplay();
}
void MyInit() {
    glClearColor(0, 0, 0, 1);  //background color
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-1, 1, -1, 1, 2, 1000);
    glMatrixMode(GL_MODELVIEW);
}
void update(int value) {
    sz += .2f;
   sharkX -=.1f;
    bombY -=.1f;
     _angle += 5.0f;
    if (_angle > 360) {
        _angle -= 360;
    }
    if(sharkX<=-2)
    {
        sharkX=2;
    }
    if(bombY<=-5)
    {
        bombY=2;
    }
    if(bx>=1.5 && bx<=1.8)
    {
        if((int(sz-5)%15)>=-2 && (int(sz-5)%15)<=2)
    {
        if(!sc)
        {
            ct++;
            sc=true;
        }
    }

        cout<<"Score :"<<ct<<endl;
    }
    else if(bx>=-0.9 && bx<=-.6)
    {
        if((int(sz+15)%15)>=-2 && (int(sz+15)%15)<=2)
    {
        if(!sc)
        {
            ct++;
            sc=true;
        }
    }

        cout<<"Score :"<<ct<<endl;
    }
    else if(bx>=-1.8 && bx<=-1.5)
    {
        if((int(sz-5)%15)>=-2 && (int(sz-5)%15)<=2)
    {
        if(!sc)
        {
            ct++;
            sc=true;
        }
    }

        cout<<"Score :"<<ct<<endl;
    }
    else{
        sc=false;
    }
    if(bx==.00001 || bx==.9 && (int(sz-5)%10)>=-1 && (int(sz-6)%10)<=1 && bombY>=-1 && bombY<=.3){
        cout<<"Final Score :"<<ct<<endl;
        exit(1);
    }
    if(bx==sharkX && (int(sz)%30>=-5) && (int(sz)%30<=5)  )   {
        cout<<"Final Score :"<<ct<<endl<<"Collision With Shark"<<endl;;
        exit(1);
    }
    if(bx>=-.9 && bx<=-.6)
    {
        if((int(sz+2)%20)>=-5 && (int(sz+2)%20)<=5)
    {
        cout<<"Final Score :"<<ct<<endl<<"Collision With Boat"<<endl;;
        exit(1);
    }
    }
    else if(bx>=.9 && bx<=1.2)
    {
          if((int(sz)%30)>=-5 && (int(sz)%30)<=5)
    {
        cout<<"Final Point :"<<ct<<endl<<"Collision With Boat"<<endl;
        exit(1);
    }
    }






    glutPostRedisplay();  ////Tell GLUT that the scene has changed
    glutTimerFunc(25, update, 0);
}
void LoadTexture(const char*filename,GLuint num)
{
    glGenTextures(1, &ID);
    glBindTexture(GL_TEXTURE_2D, ID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowPosition(50, 50);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Pirates of the Bengal");
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\download.bmp",1);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\water2.bmp",2);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\images.bmp",3);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\shark.bmp",4);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\tree.bmp",5);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\gift3.bmp",6);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\sky_front.bmp",7);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\tree_wood.bmp",8);
    LoadTexture("C:\\Users\\Dipto\\Desktop\\boat2\\bomb2.bmp",9);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    MyInit();
    glutDisplayFunc(display);
    glutKeyboardFunc(myKeyboardFunc);
    //glutReshapeFunc(resize);
    glEnable(GL_LIGHTING);
    glutTimerFunc(25, update, 0);  //Add a timer
    printf("Follow below instructions:\n\n x/X to move camera on X axis \n y/Y to move camera on Y axis \n R to look right side \n L to look left side \n U to look up \n D to look down\n n for boat moving forward\n m for boat moving backward\n b for moving boat in right direction\n v for moving boat in left direction\n 1 for on/off Right light,2 for on/off left light,3 for on/off spot light \n a for on/off ambient light \n s for on/off specular light \n d for on/off diffusion light \n");
    glutMainLoop();
    return 0;
}

