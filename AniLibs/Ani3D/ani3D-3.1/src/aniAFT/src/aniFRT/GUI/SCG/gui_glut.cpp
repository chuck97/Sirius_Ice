#include <GL/freeglut.h>
#include <GL/glext.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "gl2ps.h"
#include "mesh.h"

const char *version =
"gui_glut - 3D surface mesh generator (glut version)\n"
"\n"
"Usage: gui_glut.exe [+csp] [params]"
"\n"
"Press 'h' for help\n\n";

const char *help =
"\"+s r x y z r p q M\" - M-mesh for a sphere with radius r, centered at (x, y, z) and rotated by angles r, p, q\n"
"\"+c r h x y z r p q M\" - M-mesh for a cylinder with radius r, height h, centered at (x, y, z) and "
"rotated by angles r, p, q\n"
"\"+p a b c x y z r p q M\" - M-mesh for a parallelepiped with radius r, centered at (x, y, z) and "
"rotated by angles r, p, q\n\n"

"Usage example:\n"
"./gui_glut.exe +s 1.6 0 0 1 0 0 0 0.3 +c 0.6 5 0 0 -1 0.3 0 0.5 0.3\n"
"   creates a 0.3-mesh for sphere with radius 1.6, centered at (0, 0, 1) and rotated by angles 0, 0, 0\n"
"   and a 0.3-mesh for cylinder with radius 0.6, height 5, centered at (0, 0, -1) and rotated by angles 0.3, 0 and 0.5\n"
"\n"
"./gui_glut.exe +p 1.3 1.5 0.6 0 0 1 0.2 0 0 0.2 +c 0.6 5 0 0 -1 0.3 0 0.5 0.2\n"
"   creates a 0.2-mesh for parallelepiped with sides 1.3, 1.5 and 0.6, centered at (0, 0, 1) and \n"
"   rotated by angles 0.2, 0 and 0\n"
"   and a 0.2-mesh for cylinder with radius 0.6, height 5, centered at (0, 0, -1) and rotated by angles 0.3, 0 and 0.5\n"

"W:w write info frt-gui.frt, then you can use ./main_aft.exe frt-gui.frt to build tetra\n"
"Esc,q:\texit\n"
"h:\tthis help\n"
"left  mouse button, arrow keys: rotate\n"
"right mouse button, w,a,s,d: move\n"
"Z,z: zoom\n"
"n:\tunite\n"
"N:\tintersect\n"
"M,m:\tsubstract\n"
"V:\tdraw points\n"
"v:\tdraw vertices\n"
"E:\tdraw edges\n"
"e:\tdraw face edges\n"
"f:\tdraw faces\n"
"o:\tback offset\n"
"I,i:\tline width\n"
"U,u:\tpoint size\n"
"P:\tsave screenshot to eps file\n"
"p:\tsave screenshot to pdf file\n"
"C:\tsave camera position and direction\n"
"c:\tload camera position and direction\n"   
"B:\tbenchmark\n";

/* lighting parameters */
double back_color[4] = {1.0, 1.0, 1.0, 0.0};
double face_color[4] = {0.5, 0.5, 0.5, 0.9};
double edge_color[4] = {0.0, 0.0, 0.0, 0.9};
double vert_color[4] = {0.8, 0.8, 0.8, 0.9};
float light_ambt[4] = {0.5, 0.5, 0.5, 1.0};
float light_diff[4] = {0.9, 0.9, 0.9, 1.0};

float light_pos[4] = {-2, 1, 3, 0};

/* geometry to be rendered */
typedef double vertex[3];
typedef unsigned tri[3];
typedef unsigned edge[2];
typedef unsigned point;
int     n_vertices;
int     n_tris;
int     n_edges;
tri    *tris;
edge   *edges;
vertex *vertices;
vertex *normals;

mesh** msh;
int mesh_max_num = 10;

int current = -1, num_total = -1;


/* renderer vars */
double cur_mv_matrix[4][4];
double cur_zoom = 2.0;
int    cur_vp_x, cur_vp_y;

/* state vars */
int benchmarking     = 0;
int draw_vertices    = 0;
int draw_edges       = 1;
int draw_faces       = 1;
int draw_faceedges   = 1;
int light_screen     = 1;
int back_offset      = 0;
int both_sides       = 0;

double LineWidth = 1.0;
double PointSize = 5.0;
double PaperMult = 1.0;

/* performance vars */
unsigned frames_displayed;
int time_start;

/* filename for output ps/pdf */
char fname[1024];

#define checkGetProcAddress(proc) assert(proc = glutGetProcAddress(#proc));


static void normalize(GLdouble *n) 
{
	double s = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	if (s != 0) s = 1/s;
	n[0] *= s, n[1] *= s, n[2] *= s;
}

static void mknorm(GLdouble *a, GLdouble *b, GLdouble *c, GLdouble *n) 
{
	vertex u, v;
	u[0] = b[0]-a[0]; v[0] = c[0]-a[0];
	u[1] = b[1]-a[1]; v[1] = c[1]-a[1];
	u[2] = b[2]-a[2]; v[2] = c[2]-a[2];
	n[0] = u[1]*v[2] - v[1]*u[2];
	n[1] = u[2]*v[0] - v[2]*u[0];
	n[2] = u[0]*v[1] - v[0]*u[1];
	normalize(n);
}


static void change_lighting(void) 
{
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diff);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_ambt);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, both_sides);
	glClearColor(back_color[0], back_color[1], back_color[2], back_color[3]);
}

static void setup_lighting(void) 
{
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	change_lighting();
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glEnable(GL_COLOR_MATERIAL);
}

static void setup_offset(void) 
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	if (back_offset) 
		glPolygonOffset(-1.0, -1.0);
	else 
		glPolygonOffset(1.0, 1.0);
	gl2psEnable(GL2PS_POLYGON_OFFSET_FILL);
}

static void setup_depth(void) 
{
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
}

static void setup_size(void)
{
	glLineWidth(LineWidth);
	glPointSize(PointSize);
}

static int makeObject (int argc, char *argv[])
{
    int i, k, ret;;
    double param[10];
    msh = (mesh**)malloc (sizeof (mesh*) * mesh_max_num);
    for (i = 0; i < mesh_max_num; i++)
        msh[i] = new mesh ();

    if (argc == 1)
    {
	fprintf(stderr, "%s", help);
	exit(0);
    }
    for (i = 1; i < argc; i++)
    {
	if (argv[i][0]=='+') 
	{
	    switch (argv[i][1]) 
	    {
		case 'c': 
		    for (k = i + 1; k < i + 10 && k < argc; k++)
		    {
			param[k - (i + 1)] = (double)atof (argv[k]);
//			printf ("%lf\n", param[k - (i + 1)]);
		    }
		    num_total++;
		    if ((ret = msh[num_total]->primitive_cylinder (param[0], param[1], param[8])) < 0)
		    {
			printf ("Oops, %d!\n", ret);
			return 0;
		    }
		    msh[num_total]->move (param[2], param[3], param[4], param[5], param[6], param[7]);
		    i = k - 1;
		    continue;
			
		case 's':  
		    for (k = i + 1; k < i + 9 && k < argc; k++)
		    {
			param[k - (i + 1)] = (double)atof (argv[k]);
//			printf ("%lf\n", param[k - (i + 1)]);
		    }
		    num_total++;
		    if ((ret = msh[num_total]->primitive_sphere (param[0], param[7])) < 0)
		    {
			printf ("Oops, %d!\n", ret);
			return 0;
		    }
		    msh[num_total]->move (param[1], param[2], param[3], param[4], param[5], param[6]);
		    i = k - 1;
		    continue;

		case 'p': 
		    for (k = i + 1; k < i + 11 && k < argc; k++)
		    {
			param[k - (i + 1)] = (double)atof (argv[k]);
//			printf ("%lf\n", param[k - (i + 1)]);
		    }
		    num_total++;
		    if ((ret = msh[num_total]->primitive_paral (param[0], param[1], param[2], param[9])) < 0)
		    {
			printf ("Oops, %d!\n", ret);
			return 0;
		    }
		    msh[num_total]->move (param[3], param[4], param[5], param[6], param[7], param[8]);
		    i = k - 1;
		    continue;

		default:
		    continue;
	    }
	} 
    }

/*  int ret;
    if ((ret = msh[num_total]->primitive_sphere (1.5, 0.2)) < 0)
//    if ((ret = msh[num_total]->primitive_cylinder (0.9, 2, 0.2)) < 0)
//    if ((ret = msh[num_total]->primitive_paral (1, 1.2, 1, 0.2)) < 0)
    {
        printf ("Oops, %d!\n", ret);
        return 0;
    }
    num_total++;

//    if ((ret = msh[num_total]->primitive_sphere (2, 0.2)) < 0)
    if ((ret = msh[num_total]->primitive_cylinder (0.5, 5, 0.2)) < 0)
//    if ((ret = msh[num_total]->primitive_paral (1, 1.5, 2, 0.2)) < 0)
    {
        printf ("Oops, %d!\n", ret);
        return 0;
    }
    
    msh[num_total]->move (0.6, 0.5, -2.3, 0.2, 0.1, 0);
    num_total++; */

    return 1;
}


static void init(int argc, char *argv[])
{
        (void)argc;
        (void)argv;
	setup_lighting();
	setup_depth();
	setup_offset();
	glShadeModel(GL_FLAT);
	/*Enable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);*/
	setup_size();	
	glClearColor(back_color[0], back_color[1], back_color[2], back_color[3]);
	glEnableClientState (GL_VERTEX_ARRAY);
}

void setup_ortho(void) 
{
	int w = cur_vp_x;
	int h = cur_vp_y;
	double z = cur_zoom;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w>h)
		glOrtho(-z*w/h, z*w/h, -z, z, -100, 100);
	else
		glOrtho(-z, z, -z*h/w, z*h/w, -100, 100);
	glMatrixMode(GL_MODELVIEW);
}

void reshape(int w, int h) 
{
	glViewport(0, 0, w, h);
	cur_vp_x = w;
	cur_vp_y = h;
	setup_ortho();
}

void display(void) 
{
	int i, m_num;
	GLdouble normals[3];
	glClear(GL_COLOR_BUFFER_BIT | (GL_DEPTH_BUFFER_BIT));

	if (draw_faces)
	{
	    glColor4dv(face_color);
	    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	    for (m_num = 0; m_num <= num_total; m_num++)
	    {
		if (msh[m_num]->hide)
		    continue;

		glVertexPointer (3, GL_DOUBLE, 0, msh[m_num]->Vertex);

		for (i = 0; i < msh[m_num]->nF; i++) 
		{
		    mknorm (msh[m_num]->Vertex + 3 * msh[m_num]->Index[i * 3], 
			    msh[m_num]->Vertex + 3 * msh[m_num]->Index[i * 3 + 1], 
			    msh[m_num]->Vertex + 3 * msh[m_num]->Index[i * 3 + 2], 
			    normals);
		    glNormal3dv(normals);
		    glDrawElements (GL_TRIANGLES, 3, GL_UNSIGNED_INT, msh[m_num]->Index + i * 3);
		}
	    }
	}
		
	if (draw_faceedges) 
	{
	    glColor4dv(edge_color);
	    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	    for (m_num = 0; m_num <= num_total; m_num++)
	    {
		if (msh[m_num]->hide)
		    continue;

                glVertexPointer (3, GL_DOUBLE, 0, msh[m_num]->Vertex);

		glDrawElements (GL_TRIANGLES, 3 * msh[m_num]->nF, GL_UNSIGNED_INT, msh[m_num]->Index);
	    }
	}
		
	if (draw_vertices) 
	{
	    glColor4dv(vert_color);
	    for (m_num = 1; m_num <= num_total; m_num++)
	    {
		if (msh[m_num]->hide)
		    continue;

		glVertexPointer (3, GL_DOUBLE, 0, msh[m_num]->Vertex);

		for (i = 0; i < msh[m_num]->nV; i++)	
		    glDrawElements(GL_POINTS, 1, GL_UNSIGNED_INT, &i);
	    }
	}
		
	glutSwapBuffers();
	glutReportErrors();
	frames_displayed++;
}

void idle() 
{
	glRotatef(1, 1, 0, 0);
	glutPostRedisplay();
}

void save_screenshot(int save_to_pdf) 
{
	char fn[1024];
	FILE *outfile;
	int buffsize = 0;
	static int num[2]={0, 0};
	if (num[save_to_pdf]) sprintf(fn, (save_to_pdf)?"%s-%d.pdf":"%s-%d.eps", fname, num[save_to_pdf]);
	else sprintf(fn, (save_to_pdf)?"%s.pdf":"%s.eps", fname);
	outfile = fopen(fn, "wb");
	if (!outfile) 
	{
		perror(fn);
		return;
	}
	
	if (!save_to_pdf) 
	{
		gl2psBlendFunc(GL_ONE, GL_ZERO);
		gl2psDisable(GL2PS_BLEND);
	}
	
	do 
	{
		buffsize += 8*1024*1024;
		gl2psBeginPage("Mesh", "SMV+gl2ps", 0, save_to_pdf?GL2PS_PDF:GL2PS_PS, GL2PS_BSP_SORT,
				GL2PS_NO_PS3_SHADING | /*GL2PS_DRAW_BACKGROUND | */
				GL2PS_USE_CURRENT_VIEWPORT  | 
				((save_to_pdf)?GL2PS_BEST_ROOT:GL2PS_OCCLUSION_CULL),
				GL_RGBA, 0, 0, 0, 0, 0, buffsize, outfile, 0);
		setup_offset();
		gl2psLineWidth(LineWidth*PaperMult);
		gl2psPointSize(PointSize*PaperMult);
		display();
	} 
	
	while (gl2psEndPage() == GL2PS_OVERFLOW);

	fclose(outfile);
	printf("%s\n", fn);
	num[save_to_pdf]++;
}

void keyboard_motion(unsigned char key) 
{
	double sens = 0.05*cur_zoom;
	if (key == tolower(key))
		sens *= 8;
	key = tolower(key);
	glGetDoublev(GL_MODELVIEW_MATRIX, &cur_mv_matrix[0][0]);
	glLoadIdentity();
	switch (key) 
	{
		case 'w':
			glTranslatef(0, sens, 0);
			break;
		case 's':
			glTranslatef(0, -sens, 0);
			break;
		case 'a':
			glTranslatef(-sens, 0, 0);
			break;
		case 'd':
			glTranslatef(sens, 0, 0);
			break;
	}
	glMultMatrixd(&cur_mv_matrix[0][0]);
}

void save_cam(void) 
{
	char fn[1024];
	FILE *camfile;
	int i, j;
	sprintf(fn, "%s.cam", fname);
	camfile = fopen(fn, "w");
	if (!camfile) 
	{
		perror(fn);
		return;
	}
	fprintf(camfile, "%.4lf\n\n", cur_zoom);
	glGetDoublev(GL_PROJECTION_MATRIX, &cur_mv_matrix[0][0]);
	for (i=0; i<4; i++) 
	{
		for (j=0; j<4; j++)
			fprintf(camfile, "%.4lf ", cur_mv_matrix[i][j]);
		fprintf(camfile, "\n");
	}
	fprintf(camfile, "\n");
	glGetDoublev(GL_MODELVIEW_MATRIX, &cur_mv_matrix[0][0]);
	for (i=0; i<4; i++) 
	{
		for (j=0; j<4; j++)
			fprintf(camfile, "%.4lf ", cur_mv_matrix[i][j]);
		fprintf(camfile, "\n");
	}
	fclose(camfile);
}

void load_cam(void) 
{
	char fn[1024];
	FILE *camfile;
	int i, j;
	sprintf(fn, "%s.cam", fname);
	camfile = fopen(fn, "r");
	if (!camfile) 
	{
		perror(fn);
		return;
	}
	fscanf(camfile, "%lf", &cur_zoom);
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			fscanf(camfile, "%lf", &cur_mv_matrix[i][j]);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&cur_mv_matrix[0][0]);

	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			fscanf(camfile, "%lf ", &cur_mv_matrix[i][j]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&cur_mv_matrix[0][0]);
	fclose(camfile);
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 27:
		case 'q':
		case 'Q':
			exit(0);
		case 'h':
			fprintf(stderr, "%s", help);
			break;
		case 'B':
			if (benchmarking) {
				printf("%g\n", frames_displayed*1000.0/(glutGet(GLUT_ELAPSED_TIME)-time_start));
				glutIdleFunc(0);
			} else {
				frames_displayed = 0;
				time_start = glutGet(GLUT_ELAPSED_TIME);
				glutIdleFunc(idle);
			}
			benchmarking ^= 1;
			break;
		case 'o':
			back_offset ^= 1;
			setup_offset();
			break;
		case 'E':
			draw_edges ^= 1;
			break;
		case 'e':
			draw_faceedges ^= 1;
			break;
		case 'f':
			draw_faces ^= 1;
			break;
		case 'v':
			draw_vertices ^= 1;
			break;
		case 'I':
		case 'i':
			LineWidth *= (key=='I')?(1.25):(0.8);
			setup_size();
			break;
		case 'U':
		case 'u':
			PointSize *= (key=='U')?(1.25):(0.8);
			setup_size();
			break;
		case 'p':
			save_screenshot(1);
			break;
		case 'P':
			save_screenshot(0);
			break;
		case 'C':
			save_cam();
			break;
		case 'c':
			load_cam();
			break;
		case 'N':
			msh[0]->in = 1;
			msh[1]->in = 1;
			msh[0]->intersect (msh[1]);
			break;
		case 'n':
			msh[0]->in = 0;
			msh[1]->in = 0;
			msh[0]->intersect (msh[1]);
			break;
		case 'M':
			msh[0]->in = 1;
			msh[1]->in = 0;
			msh[0]->intersect (msh[1]);
			break;
		case 'm':
			msh[0]->in = 0;
			msh[1]->in = 1;
			msh[0]->intersect (msh[1]);
			break;
		case 'w':
		case 'W':
			msh[0]->mesh_write_smv("frt-gui.frt");
			break;
                case 'z':
                   	cur_zoom *= 1-0.05;
	                setup_ortho();
	                glutPostRedisplay();
                        break;
	        case 'Z':
                   	cur_zoom *= 1+0.05;
	                setup_ortho();
	                glutPostRedisplay();
                        break;
		case 'a':
		case 'd':
		case 's':
		case 'A':
		case 'D':
		case 'S':


			keyboard_motion(key);
			break;
	}
	glutPostRedisplay();
	(void)x, (void)y;
}

void keyboard_special(int key, int x, int y) 
{
	double sens = 5;
	glGetDoublev(GL_MODELVIEW_MATRIX, &cur_mv_matrix[0][0]);
	glLoadIdentity();
	switch (key) 
	{
		case GLUT_KEY_LEFT:
			glRotatef(-sens, 0, 1, 0);
			break;
		case GLUT_KEY_RIGHT:
			glRotatef(sens, 0, 1, 0);
			break;
		case GLUT_KEY_DOWN:
			glRotatef(sens, 1, 0, 0);
			break;
		case GLUT_KEY_UP:
			glRotatef(-sens, 1, 0, 0);
			break;
	}
	glMultMatrixd(&cur_mv_matrix[0][0]);
	glutPostRedisplay();
	(void) x; (void) y;
}

int motion_start_x, motion_start_y, m_button;

void mouse_key(int button, int state, int x, int y) 
{
	if (state == GLUT_DOWN) 
	{
		motion_start_x = x;
		motion_start_y = y;
		glGetDoublev(GL_MODELVIEW_MATRIX, &cur_mv_matrix[0][0]);
		m_button = button;
		if (button==3 || button==4) {
		    cur_zoom *= 1+((button==3) ? 1 : -1)*0.05;
		    setup_ortho();
		    glutPostRedisplay();
		}
	}
}

void mouse_motion(int x, int y) 
{
	glLoadIdentity();
	switch (m_button) 
	{
		case GLUT_LEFT_BUTTON:
			glRotatef(y-motion_start_y, 1, 0, 0);
			glRotatef(x-motion_start_x, 0, 1, 0);
			break;
		case GLUT_RIGHT_BUTTON:
			glTranslatef((x-motion_start_x)*0.01*cur_zoom, (motion_start_y-y)*0.01*cur_zoom, 0);
			break;
	}
	glMultMatrixd(&cur_mv_matrix[0][0]);
	glutPostRedisplay();
}

void mouse_wheel(int button, int state, int x, int y) 
{
	cur_zoom *= 1+state*0.05;
	setup_ortho();
	glutPostRedisplay();
	(void) button; (void) x; (void) y;
}

int main(int argc, char *argv[]) 
{
	int param=1, i, j, needcam=0, makeshot=0;
	char *c;
	for (i=1; i<argc; i++)
	{
		if ((param) && (argv[i][0]=='-')) 
		{
			for (j=1; argv[i][j]; j++)
				switch (argv[i][j]) 
				{
					case 'v': draw_vertices ^= 1; break;
					case 'E': draw_edges ^= 1; break;
					case 'e': draw_faceedges ^= 1; break;
					case 'f': draw_faces ^= 1; break;
					case 'k': both_sides ^= 1; break;
					case 'l': light_screen ^= 1; break;
					case 'o': back_offset ^= 1; break;
					case 'U': PointSize *= 1.25; break;
					case 'u': PointSize *= 0.8; break;
					case 'I': LineWidth *= 1.25; break;
					case 'i': LineWidth *= 0.8; break;
					case 'P': makeshot ^= 1; break;
					case 'p': makeshot ^= 2; break;
					case 'c': needcam ^= 1; break;
					case '-': param=0; break;
				}
		} 
		else break;
	}
		
	argc -= i-1; argv += i-1;
	if (argc<2) strcpy(fname, "null");
	if (argc>2) strcpy(fname, "mix");
	if (argc==2) strcpy(fname, argv[1]);
	c = fname+strlen(fname)-1;
	if ((c>fname) && (*c--=='v')) if ((c>fname) && (*c--=='m'))
		if ((c>fname) && (*c--=='s')) if ((c>fname) && (*c=='.'))
			*c = 0;
	if (argc<2) printf("%s", version);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(790, 610);
	glutCreateWindow("mesh_3D");
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutSpecialUpFunc(keyboard_special);
	glutMotionFunc(mouse_motion);
	glutMouseFunc(mouse_key);
	glutMouseWheelFunc(mouse_wheel);
	init(argc, argv);
	makeObject(argc, argv);
	if (needcam) load_cam();
	if (makeshot) 
	{
		if (makeshot & 1) save_screenshot(0);
		if (makeshot & 2) save_screenshot(1);
		return 0;
	}
	glutMainLoop();
	return 0;
}


