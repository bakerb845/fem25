#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <visit_writer.h>

/*****************************************************************************/
#ifdef HPUX
void vtk_model(filnam,LENFL, LDD,NNPG,NELEM,NIENGV, 
               iengnod,iengv, xpts,zpts,dat)
#else
void vtk_model_(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,
                iengnod,iengv, xpts,zpts,dat)
#endif
/*
   vtk_model: Writes the isotropic Vp, Vs, and Density model for a 2D medium
*/
int *LENFL;    //Length of filename
char *filnam;  //Filename
int *LDD;      //Leading dimension of resp 
int *NELEM;    //Number of elements
int *NNPG;     //Number of anchor nodes in mesh 
int *NIENGV;   //Size of iengv
int *iengv;    //Vector form IENG connectivity 
int *iengnod;  //3 or 4 corresponding to anchor nodes for tris or quads
float *xpts;   //X anchor node locations 
float *zpts;   //Z anchor node locations
float *dat;    //Vp velocity, Vs Velocity, Density 
               //packed as [0:nnpg-1,nnpg:2*nnpg-1:2*nnpg:3*nnpg-1]
{
   // Variables for model mesh 
   int nnodes = *NNPG; int nnodes3 = 3*nnodes, ldd = (int)*LDD;
   float pts[nnodes3];        //Holds element points 
   int nzones = *NELEM;       //Number of elements
   int zonetypes[nzones];     //Tells VisIt if element is a quad or triangle
   int nsc = *NIENGV;         //Size of connectivity
   int connectivity[nsc];     //Vector form of IENG with C numbering
   float vp[nnodes];          //Compressional velocity 
   float vs[nnodes];          //Shear velocity
   float rho[nnodes];         //Density
   int nvars = 3;             //Vp, Vs, Density 
   int centering[] = {1,1,1}; //No centering, info is node based 
   const char *varnames[] = {"Vp_m/s","Vs_m/s","Dens_kg/m**3"};
   int vardims[] = {1,1,1};  //Scalars 
   float *vars[] = {vp,vs,rho}; //Vp vs and rho 
   // Variables for filename
   int lenfl = *LENFL;
   char filename[lenfl+1];
   // Set points array and data 
   int ipts = 0; int inpg, loc1, loc2, loc3;
   for (inpg=0; inpg<nnodes; ++inpg){
       pts[ipts+0] = xpts[inpg];
       pts[ipts+1] = 0.; 
       pts[ipts+2] = zpts[inpg];
       loc1 = inpg; 
       loc2 = ldd + inpg;
       loc3 = 2*ldd + inpg;
       vp[inpg] = dat[loc1];
       vs[inpg] = dat[loc2];
       rho[inpg] = dat[loc3]; 
       ipts = ipts + 3;
   }   
   // Set zone types and connectivity
   int isc = 0;
   int ielem, nloop, ia;
   for (ielem=0; ielem<nzones; ++ielem){
       if (iengnod[ielem] == 4){
          nloop = 4;
          zonetypes[ielem] = VISIT_QUAD;
       }
       else if (iengnod[ielem] == 3){
          nloop = 3;
          zonetypes[ielem] = VISIT_TRIANGLE;
       }
       else{
          printf(" vtk_model: Error element type undefined!\n");
          return;
       }
       for (ia=0; ia<nloop; ++ia){
           connectivity[isc] = iengv[isc] - 1; //Fortran -> C
           isc = isc + 1;
       }
   }
   // Set filename
   strncpy(filename,filnam,lenfl);
   filename[lenfl] = '\0';
   printf(" vtk_model: Writing %s \n",filename);
   write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                     zonetypes, connectivity, nvars, vardims, centering,
                     varnames, vars);
   bzero(filename,lenfl);
   return;
}
/*****************************************************************************/
#ifdef HPUX
void vtk_shgrad2(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,ITYPE, 
                 iengnod,iengv, xpts,zpts,dat)
#else
void vtk_shgrad2_(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,ITYPE,
                  iengnod,iengv, xpts,zpts,dat)
#endif
/*
   vtk_shgrad2: Writes the Vp and Vs search direction, hessian, or gradient
*/
int *LENFL;    //Length of filename
char *filnam;  //Filename
int *LDD;      //Leading dimension of resp 
int *NELEM;    //Number of elements
int *NNPG;     //Number of anchor nodes in mesh 
int *NIENGV;   //Size of iengv
int *ITYPE;    //1 -> Gradient, 2 -> Hessian, 3 -> Search
int *iengv;    //Vector form IENG connectivity 
int *iengnod;  //3 or 4 corresponding to anchor nodes for tris or quads
float *xpts;   //X anchor node locations 
float *zpts;   //Z anchor node locations
float *dat;    //Gp velocity and Gs Velocity or Hp velocity and Hs velocity  or
               //Ss velocity and Ss velocity search direction
               //packed as [0:nnpg-1,nnpg:2*nnpg-1]
{
   // Variables for model mesh 
   int nnodes = *NNPG; int nnodes3 = 3*nnodes, ldd = (int)*LDD;
   float pts[nnodes3];        //Holds element points 
   int nzones = *NELEM;       //Number of elements
   int zonetypes[nzones];     //Tells VisIt if element is a quad or triangle
   int nsc = *NIENGV;         //Size of connectivity
   int connectivity[nsc];     //Vector form of IENG with C numbering
   float hgp[nnodes];         //P gradient or P hessian
   float hgs[nnodes];         //S gradient or S hessian
   int nvars = 2;             //Vp, Vs, Density 
   int centering[] = {1,1}; //No centering, info is node based 
   const char *varnamesG[] = {"GradP_m/s","Grads_m/s"};
   const char *varnamesH[] = {"HessianP","HessianS"};
   const char *varnamesS[] = {"dP_m/s","dS_m/s"};
   int vardims[] = {1,1};  //Scalars 
   float *vars[] = {hgp,hgs}; //Vp vs and rho 
   // Variables for filename
   int lenfl = *LENFL;
   char filename[lenfl+1];
   // Set points array and data 
   int ipts = 0; int inpg, loc1, loc2;
   int ivinv;
   for (inpg=0; inpg<nnodes; ++inpg){
       pts[ipts+0] = xpts[inpg];
       pts[ipts+1] = 0.;
       pts[ipts+2] = zpts[inpg];
       loc1 = inpg;
       loc2 = nnodes + inpg;
       hgp[inpg] = dat[loc1];
       hgs[inpg] = dat[loc2];
       ipts = ipts + 3;
   }
   // Set zone types and connectivity
   int isc = 0;
   int ielem, nloop, ia;
   for (ielem=0; ielem<nzones; ++ielem){
       if (iengnod[ielem] == 4){
          nloop = 4;
          zonetypes[ielem] = VISIT_QUAD;
       }
       else if (iengnod[ielem] == 3){
          nloop = 3;
          zonetypes[ielem] = VISIT_TRIANGLE;
       }
       else{
          printf(" vtk_hgrad2: Error element type undefined!\n");
          return;
       }
       for (ia=0; ia<nloop; ++ia){
           connectivity[isc] = iengv[isc] - 1; //Fortran -> C
           isc = isc + 1;
       }
   }
   // Set filename
   strncpy(filename,filnam,lenfl);
   filename[lenfl] = '\0';
   printf(" vtk_hgrad2: Writing %s \n",filename);
   if (*ITYPE == 1){
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesG, vars);
   }
   else if (*ITYPE == 2) {
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesH, vars);
   }
   else if (*ITYPE == 3) {
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesS, vars);
   }  
   else {
      printf(" vtk_shgrad2: I dont know what to plot\n");
   }
   bzero(filename,lenfl);
   return;
}
/*****************************************************************************/
#ifdef HPUX
void vtk_shgrad1(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,ITYPE, 
                 iengnod,iengv, xpts,zpts,dat)
#else
void vtk_shgrad1_(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,ITYPE,
                  iengnod,iengv, xpts,zpts,dat)
#endif
/*
   vtk_shgrad1: Writes Vp or Vs search direction, hessian, or gradient
*/
int *LENFL;    //Length of filename
char *filnam;  //Filename
int *LDD;      //Leading dimension of resp 
int *NELEM;    //Number of elements
int *NNPG;     //Number of anchor nodes in mesh 
int *NIENGV;   //Size of iengv
int *ITYPE;    //1 -> P Gradient, 2 -> S gradient, 
               //3 -> P Hessian,  4 -> S Hessian, 
               //5 -> P search,   6 -> S search direction 
int *iengv;    //Vector form IENG connectivity 
int *iengnod;  //3 or 4 corresponding to anchor nodes for tris or quads
float *xpts;   //X anchor node locations 
float *zpts;   //Z anchor node locations
float *dat;    //Gp velocity, Gs Velocity, Hp hessian or Hs hessian
{
   // Variables for model mesh 
   int nnodes = *NNPG; int nnodes3 = 3*nnodes, ldd = (int)*LDD;
   float pts[nnodes3];        //Holds element points 
   int nzones = *NELEM;       //Number of elements
   int zonetypes[nzones];     //Tells VisIt if element is a quad or triangle
   int nsc = *NIENGV;         //Size of connectivity
   int connectivity[nsc];     //Vector form of IENG with C numbering
   float gh[nnodes];          //Gradeint or Hessian 
   int nvars = 1;             //Vp, Vs, Density 
   int centering[] = {1}; //No centering, info is node based 
   const char *varnamesGP[] = {"GradP_m/s"};
   const char *varnamesGS[] = {"GradS_m/s"};
   const char *varnamesHP[] = {"HessianP"};
   const char *varnamesHS[] = {"HessianS"};
   const char *varnamesHSP[]= {"HessianPS"};
   const char *varnamesSP[] = {"dP_m/s"};
   const char *varnamesSS[] = {"dS_m/s"};

   int vardims[] = {1};  //Scalars 
   float *vars[] = {gh}; //Gradient or Hessian
   // Variables for filename
   int lenfl = *LENFL;
   char filename[lenfl+1];
   // Set points array and data 
   int ipts = 0; int inpg, loc1;
   int ivinv;
   for (inpg=0; inpg<nnodes; ++inpg){
       pts[ipts+0] = xpts[inpg];
       pts[ipts+1] = 0.;
       pts[ipts+2] = zpts[inpg];
       loc1 = inpg;
       gh[inpg] = dat[loc1];
       ipts = ipts + 3;
   }
   // Set zone types and connectivity
   int isc = 0;
   int ielem, nloop, ia;
   for (ielem=0; ielem<nzones; ++ielem){
       if (iengnod[ielem] == 4){
          nloop = 4;
          zonetypes[ielem] = VISIT_QUAD;
       }
       else if (iengnod[ielem] == 3){
          nloop = 3;
          zonetypes[ielem] = VISIT_TRIANGLE;
       }
       else{
          printf(" vtk_shgrad1: Error element type undefined!\n");
          return;
       }
       for (ia=0; ia<nloop; ++ia){
           connectivity[isc] = iengv[isc] - 1; //Fortran -> C
           isc = isc + 1;
       }
   }
   // Set filename
   strncpy(filename,filnam,lenfl);
   filename[lenfl] = '\0';
   printf(" vtk_shgrad1: Writing %s \n",filename);
   if (*ITYPE == 1){
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesGP, vars);
   }
   else if (*ITYPE == 2){
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesGS, vars);
   }
   else if (*ITYPE == 3){
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesHP, vars);
   }
   else if (*ITYPE == 4){
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesHS, vars);
   }
   else if (*ITYPE == 5){ 
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesSP, vars);
   }   
   else if (*ITYPE == 6){ 
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesSS, vars);
   }
   else if (*ITYPE == 7){
      write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                        zonetypes, connectivity, nvars, vardims, centering,
                        varnamesHSP, vars);
   }   
   else{
      printf(" vtk_shgrad1: I dont know what to write!\n");
   }
   bzero(filename,lenfl);
   return;
}


/*****************************************************************************/
#ifdef HPUX
void vtk_rresp25d_el(filnam,LENFL, LDD,NNPG,NELEM,NIENGV, 
                     iengnod,iengv, resp,xpts,zpts)
#else
void vtk_rresp25d_el_(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,
                      iengnod,iengv, resp,xpts,zpts)
#endif
/*
   vtk_rresp25d_el: Writes real the elastic response for a 2.5D medium 
*/
int *LENFL;    //Length of filename
char *filnam;  //Filename
int *LDD;      //Leading dimension of resp 
int *NELEM;    //Number of elements
int *NNPG;     //Number of anchor nodes in mesh 
int *NIENGV;   //Size of iengv
int *iengv;    //Vector form IENG connectivity 
int *iengnod;  //3 or 4 corresponding to anchor nodes for tris or quads
float *resp;   //(u,v,w) real response at anchor nodes
float *xpts;   //X anchor node locations 
float *zpts;   //Z anchor node locations
{
   // Variables for response mesh 
   int nnodes = *NNPG; int nnodes3 = 3*nnodes, ldd = (int)*LDD;
   float pts[nnodes3];       //Holds element points 
   int nzones = *NELEM;      //Number of elements
   int zonetypes[nzones];    //Tells VisIt if element is a quad or triangle
   int nsc = *NIENGV;        //Size of connectivity
   int connectivity[nsc];    //Vector form of IENG with C numbering
   float uwresp[nnodes][3];  //Will hold horizontal response and vertical responses
   int nvars = 1;            //Only 2 variables 
   int centering[] = {1};    //No centering, info is node based 
   const char *varnames[] = {"Resp_m"};
   int vardims[] = {3};      //Vector 
   float *vars[] = {(float*)uwresp}; //Real horizontal and vertical displacements  
   // Variables for filename
   int lenfl = (int)*LENFL;
   char filename[lenfl+1];
   // Set points array and data 
   int ipts = 0; int inpg;
   for (inpg=0; inpg<nnodes; ++inpg){
       pts[ipts+0] = xpts[inpg];
       pts[ipts+1] = 0.; 
       pts[ipts+2] = zpts[inpg];
       uwresp[inpg][0] = resp[inpg];
       uwresp[inpg][1] = resp[ldd+inpg]; 
       uwresp[inpg][2] = resp[2*ldd+inpg];
       ipts = ipts + 3;
   }
   // Set zone types and connectivity
   int isc = 0; 
   int ielem, nloop, ia;
   for (ielem=0; ielem<nzones; ++ielem){
       if (iengnod[ielem] == 4){
          nloop = 4;
          zonetypes[ielem] = VISIT_QUAD;
       }
       else if (iengnod[ielem] == 3){
          nloop = 3;
          zonetypes[ielem] = VISIT_TRIANGLE;
       }
       else{
          printf(" vtk_rresp25d_el: Error element type undefined!\n");
          return;
       }
       for (ia=0; ia<nloop; ++ia){
           connectivity[isc] = iengv[isc] - 1; //Fortran -> C
           isc = isc + 1;
       }
   }
   // Set filename
   strncpy(filename,filnam,lenfl);
   filename[lenfl] = '\0';
   printf(" vtk_rresp25d_el: Writing %s \n",filename);
   write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                     zonetypes, connectivity, nvars, vardims, centering,
                     varnames, vars);
   bzero(filename,lenfl);
   return;
}
/*****************************************************************************/
#ifdef HPUX
void vtk_meshp2d(filnam,LENFL, NNPG,NELEM,NIENGV,
                    iengnod,iengv, part,xpts,zpts)
#else
void vtk_meshp2d_(filnam,LENFL, NNPG,NELEM,NIENGV,
                     iengnod,iengv, part,xpts,zpts)
#endif
/*
   vtk_meshp2d: Writes the 2d mesh partition 
*/
int *LENFL;    //Length of filename
char *filnam;  //Filename
int *NELEM;    //Number of elements
int *NNPG;     //Number of anchor nodes in mesh 
int *NIENGV;   //Size of iengv
int *iengv;    //Vector form IENG connectivity 
int *iengnod;  //3 or 4 corresponding to anchor nodes for tris or quads
float *part;   //Holds elements average partition ID
float *xpts;   //X anchor node locations 
float *zpts;   //Z anchor node locations
{
   // Variables for partition
   const char *varnames[] = {"Partition"};
   int nnodes = *NNPG; int nnodes3 = 3*nnodes;
   float pts[nnodes3];      //Holds element points 
   int nzones = *NELEM;     //Number of elements
   int zonetypes[nzones];   //Tells VisIt if element is a quad or triangle
   int nsc = *NIENGV;       //Size of connectivity
   int connectivity[nsc];   //Vector form of IENG with C numbering
   float part1[nzones];     //Will hold partition
   int nvars = 1;           //Just 1 variables in vars 
   int centering[] = {0};   //Use centering to make info element based 
   int vardims[] = {1};     //Scalars 
   float *vars[] = {part1}; //Paritition number of element at it's nodes 
   // Variables for filename
   int lenfl = (int)*LENFL;
   char filename[lenfl+1];
   // Set points array and data 
   int ipts = 0; int inpg; 
   for (inpg=0; inpg<nnodes; ++inpg){
       pts[ipts+0] = xpts[inpg];
       pts[ipts+1] = 0.; 
       pts[ipts+2] = zpts[inpg];
       ipts = ipts + 3;
   }   
   // Set zone types and connectivity
   int isc = 0; int ielem, nloop, ia;
   for (ielem=0; ielem<nzones; ++ielem){
       part1[ielem] = part[ielem];
       if (iengnod[ielem] == 4){
          nloop = 4;
          zonetypes[ielem] = VISIT_QUAD;
       }
       else if (iengnod[ielem] == 3){
          nloop = 3;
          zonetypes[ielem] = VISIT_TRIANGLE;
       }
       else{
          printf(" vtk_mesh2d: Error element type undefined!\n");
          return;
       }
       for (ia=0; ia<nloop; ++ia){
           connectivity[isc] = iengv[isc] - 1; //Fortran -> C
           isc = isc + 1;
       }
   }
   // Set filename
   strncpy(filename,filnam,lenfl);
   filename[lenfl] = '\0';
   printf(" vtk_mesh2d: Writing %s\n",filename);
   //filename[lenfl] = '\0' //May be necessary on some computers
   write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                     zonetypes, connectivity, nvars, vardims, centering,
                     varnames, vars);
   bzero(filename,lenfl);
   return;
}
/*****************************************************************************/
#ifdef HPUX
void vtk_qresp25d_el(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,
                    iengnod,iengv, xpts,zpts,uvwmag,uvwphase)
#else
void vtk_qresp25d_el_(filnam,LENFL, LDD,NNPG,NELEM,NIENGV,
                     iengnod,iengv, xpts,zpts,uvwmag,uvwphase)
#endif
/*
   vtk_qresp2d_el: Writes the complex elastic response for a 2.5D medium 
*/
int *LENFL;      //Length of filename
char *filnam;    //Filename
int *LDD;        //Leading dimension of resp 
int *NELEM;      //Number of elements
int *NNPG;       //Number of anchor nodes in mesh 
int *NIENGV;     //Size of iengv
int *iengv;      //Vector form IENG connectivity 
int *iengnod;    //3 or 4 corresponding to anchor nodes for tris or quads
float *xpts;     //X anchor node locations 
float *zpts;     //Z anchor node locations
float *uvwmag;   //Magnitude u [0:nnpg-1], v [nnpg:2*nnpg-1]; w [2*nnpg:3*nnpg-1] 
float *uvwphase; //Phase u [0:nnpg-1], v [nnpg:2*nnpg-1]; w [2*nnpg:3*nnpg-1]
{
   // Variables for partition
   int nnodes = *NNPG; int nnodes3 = 3*nnodes, ldd = (int)*LDD;
   float pts[nnodes3];       //Holds element points 
   int nzones = *NELEM;      //Number of elements
   int zonetypes[nzones];    //Tells VisIt if element is a quad or triangle
   int nsc = *NIENGV;        //Size of connectivity
   int connectivity[nsc];    //Vector form of IENG with C numbering
   float rmag [nnodes][3];   //Will hold horizontal response
   float phase[nnodes][3];   //Will hold vertical response
   int nvars = 2;            //Horizontal and vertical 
   int centering[] = {1,1};  //No centering, info is node based 
   const char *varnames[] = {"Magnitude","Phase"};
   int vardims[] = {3,3};    //Vectors 
   float *vars[] = {(float*)rmag,(float*)phase}; //Magntitude and phase 
   // Variables for filename
   int lenfl = (int)*LENFL;
   char filename[lenfl+1];
   // Set points array and data 
   int ipts = 0; int inpg;
   for (inpg=0; inpg<nnodes; ++inpg){
       pts[ipts+0] = xpts[inpg];
       pts[ipts+1] = 0.;
       pts[ipts+2] = zpts[inpg];
       rmag [inpg][0] = uvwmag[inpg];
       rmag [inpg][1] = uvwmag[ldd+inpg]; 
       rmag [inpg][2] = uvwmag[2*ldd+inpg];
       phase[inpg][0] = uvwphase[inpg];
       phase[inpg][1] = uvwphase[ldd+inpg]; 
       phase[inpg][2] = uvwphase[2*ldd+inpg];
       ipts = ipts + 3;
   }
   // Set zone types and connectivity
   int isc = 0;
   int ielem, nloop, ia;
   for (ielem=0; ielem<nzones; ++ielem){
       if (iengnod[ielem] == 4){
          nloop = 4;
          zonetypes[ielem] = VISIT_QUAD;
       }
       else if (iengnod[ielem] == 3){
          nloop = 3;
          zonetypes[ielem] = VISIT_TRIANGLE;
       }
       else{
          printf(" vtk_qresp25d_el: Error element type undefined!\n");
          return;
       }
       for (ia=0; ia<nloop; ++ia){
           connectivity[isc] = iengv[isc] - 1; //Fortran -> C
           isc = isc + 1;
       }
   }
   // Set filename
   strncpy(filename,filnam,lenfl);
   filename[lenfl] = '\0';
   printf(" vtk_qresp25d_el: Writing %s \n",filename);
   write_unstructured_mesh(filename,1,nnodes,pts, nzones,
                     zonetypes, connectivity, nvars, vardims, centering,
                     varnames, vars);
   bzero(filename,lenfl);
   return;
}
/*****************************************************************************/
#ifdef HPUX
void vtk_plot_reg_hess(filnam,LENFL, N, buffer) 
#else
void vtk_plot_reg_hess_(filnam,LENFL, N, buffer)
#endif
char *filnam;    //Filename to write
int *LENFL;      //Length of the filename
int *N;          //Number of rows/columns of Hessian matrix
float *buffer;   //Holds Hessian matrix

{
   const char *varnames[] = {"Re{adj(J)J})"};
   int n = (int)*N;
   int nwork = n*n;
   int dims[] = {n,1,n};  //Square matrix
   int nvars = {1};       //Only one variable
   //float nodal[nwork];
   int vardims[] = {1};   //These are just magnitudes 
   int centering[] = {1}; //Variables on nodal points
   float *vars[] = {buffer}; 
   // Variables for filename
   int lenfl = (int)*LENFL;
   char filename[lenfl+1];
   // copy buffer
   int ix; int iz;
   int iter = 0;
   /*
   for (ix=0; ix<n; ++ix){
      for (iz=0; iz<n; ++iz){
          nodal[iter] = buffer[iter];
          iter = iter + 1; 
      }
   }
   */
   // Set filename
   strncpy(filename,filnam,lenfl);
   filename[lenfl] = '\0';
   printf(" vtk_qresp2d_el: Writing %s \n",filename);
   
   write_regular_mesh(filename,1,dims,nvars, vardims, 
                      centering,varnames,vars);
   return;
}
