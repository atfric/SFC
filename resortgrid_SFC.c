/* Code from:
 
 Rakowsky N. & Fuchs A. Efficient local resorting techniques with space filling curves applied
 to the tsunami simulation model TsunAWI. In IMUM 2011 - The 10th International Workshop
 on Multiscale (Un-)structured Mesh Numerical Modelling for coastal, shelf and global ocean 
 dynamics, August 2011. http://hdl.handle.net/10013/epic.39576.d001
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

enum sfc {A, D, C, U};

int cmp_integer(const void *a, const void *b){

  int *i,*j;

  i = (int *)a;
  j = (int *)b;
  
  return (i[0]-j[0]);
}

int mysfc(int dim, double *lat, double *lon, int *p){
  
  int i, j, k;
  int maxlevel = 14;
  int *indx, indx_add;
  int halve;
  int *perm;
  int lc;
  double xmin, xmax, xscale, ymin, ymax, yscale, mx, ny, xhalve, yhalve;
  int nmax;
  enum sfc c;
  FILE *fptr;

  nmax = 2;
  for (i=0; i<maxlevel; i++) nmax *= 2;

  xmin = xmax = lat[0];
  ymin = ymax = lon[0];
  
  for(i = 1; i < dim; i++){
    xmin = (lat[i] < xmin) ? lat[i] : xmin;
    xmax = (lat[i] > xmax) ? lat[i] : xmax;
    ymin = (lon[i] < ymin) ? lon[i] : ymin;
    ymax = (lon[i] > ymax) ? lon[i] : ymax;
  }

  xscale = ((double) nmax)/(xmax-xmin);
  yscale = ((double) nmax)/(ymax-ymin);
  
  indx = calloc(dim, sizeof(int));
  
  for(i = 0; i < dim; i++){ 

    indx[i] = 0;

    mx = lat[i] - xmin;
    ny = lon[i] - ymin;

    xhalve = xmax-xmin;
    yhalve = ymax-ymin;
    
    c = A;
    
    indx_add = nmax*nmax;

    for(j = 0; j<maxlevel; j++){

      indx_add = indx_add / 4;

      xhalve = .5*xhalve;
      yhalve = .5*yhalve;

      /* printf("%d %d %d %d\n",halve, indx_add,m,n); */
      if(mx < xhalve)
	lc = 0;
      else{
	lc = 2;
	mx -= xhalve;
      }
      if(ny >= yhalve){
	lc++;
	ny -= yhalve;
      }      
      
      switch(c){
	
      case A:

	switch(lc){
	  
	case 0:
	  c = D;
	  break;
	case 2:
	  c = C;
	  indx[i] += 3*indx_add;
	  break;
	case 3:
	  indx[i] += 2*indx_add;
	  break;
	case 1:
	  indx[i] += indx_add;   	  
	  break;
	default:
	  printf("ERROR: falsches lc\n");	  
	  break;

	}
	break;

      case U:
	
	switch(lc){

	case 3:
	  c = C;
	  break;
	case 1:
	  c = D;
	  indx[i] += 3*indx_add;	  
	  break;
	case 0:
	  indx[i] += 2*indx_add;	  
	  break;
	case 2:
	  indx[i] += indx_add;
	  break;
	default:
	  printf("ERROR: falsches lc\n");	  
	  break;

	}
	break;
	
      case C:
	
	switch(lc){
	  
	case 3: 
	  c = U;
	  break;
	case 2:
	  c = A;
	  indx[i] += 3*indx_add;	  
	  break;
	case 0:
	  indx[i] += 2*indx_add;	  
	  break;
	case 1:
	  indx[i] += indx_add;
	  break;
	default:
	  printf("ERROR: falsches lc\n");	  
	  break;
	  
	}
	break;

      case D:
	
	switch(lc){
	  
	case 0:
	  c = A;
	  break;
	case 1:
	  c = U;
	  indx[i] += 3*indx_add;	  
	  break;
	case 3:
	  indx[i] += 2*indx_add;	  
	  break;
	case 2:
	  indx[i] += indx_add;
	  break;
	default:
	  printf("ERROR: falsches lc\n");	  
	  break;
	  
	}	  
	break;	
	  
      }	
    }    
  }

   
  perm = malloc(2*dim*sizeof(int));
  for(i = 0; i < dim; i++){
    perm[2*i] = indx[i];
    perm[2*i+1] = i;
  }  
  qsort(perm, dim, 2*sizeof(int), cmp_integer);
  
  for(i = 0; i < dim; i++)
    p[i] = perm[2*i+1];

  free(indx);
  free(perm);
}


 


int resort_grid(char *path, char *suffix){
    
  int i, j, k, l, len, a, b, c;
  int nod2D, elem2D;
  int *el, *ind, *perm, *iperm;
  double *x, *y, third = 1./3.;
  char fpath[128], fname[128];
  FILE *fptr;
    
  strcpy(fpath, path); 
  len = strlen(fpath);
  strcat(fpath, "/nod2d");
  strcat(fpath, suffix);
  
  fptr = fopen(fpath, "r");
  if(fptr == NULL){
    printf("ERROR: %s konnte nicht geoeffnet werden.\n",fpath);
    return(1);
  }
  fscanf(fptr, "%d", &nod2D);
  x    = malloc(nod2D*sizeof(double));
  y    = malloc(nod2D*sizeof(double));  
  ind  = malloc(nod2D*sizeof(int));
  for(i = 0; i < nod2D; i++)
    fscanf(fptr, "%d %lf %lf %d", &k, &x[i], &y[i], &ind[i]);  
  fclose(fptr);
  
  perm  = calloc(nod2D,sizeof(int));
  iperm = calloc(nod2D,sizeof(int));

  mysfc(nod2D, x, y, perm);
  
  for(i = 0; i < nod2D; i++)
    iperm[perm[i]] = i+1;

  strcpy(&fpath[len],"/nod2d_sfc");
  strcat(fpath, suffix);
  fptr = fopen(fpath, "w");
  fprintf(fptr, "%d\n", nod2D);
  for( i = 0; i < nod2D; i++) 
     fprintf(fptr, "%d %.15f %.15f %d\n", i+1, x[perm[i]], y[perm[i]], ind[perm[i]]);
  fclose(fptr);
    
  free(y);

  strcpy(&fpath[len], "/nodhn");
  strcat(fpath, suffix);
  
  fptr = fopen(fpath, "r");
  if(fptr != NULL){
    for(i = 0; i < nod2D; i++)
      fscanf(fptr,"%lf",&x[i]);
    fclose(fptr);
    strcpy(&fpath[len], "/nodhn_sfc");
    strcat(fpath, suffix);
    fptr = fopen(fpath, "w");
    for(i = 0; i < nod2D; i++)
      fprintf(fptr, "%lf\n", x[perm[i]]);
    fclose(fptr);
  }
  
  fptr = fopen("nodhn.out.alldata", "r");
  if(fptr != NULL){
    for(i = 0; i < nod2D; i++) fscanf(fptr,"%lf",&x[i]);
    fclose(fptr);
    fptr = fopen("nodhn_sfc.out.alldata", "w");
    for(i = 0; i < nod2D; i++) fprintf(fptr, "%lf\n", x[perm[i]]);
    fclose(fptr);
  }
 
  fptr = fopen("nodhn.out.gebco", "r");
  if(fptr != NULL){
    for(i = 0; i < nod2D; i++) fscanf(fptr,"%lf",&x[i]);
    fclose(fptr);
    fptr = fopen("nodhn_sfc.out.gebco", "w");
    for(i = 0; i < nod2D; i++) fprintf(fptr, "%lf\n", x[perm[i]]);
    fclose(fptr);
  }

  fptr = fopen("nodhn.out.orig", "r");
  if(fptr != NULL){
    for(i = 0; i < nod2D; i++) fscanf(fptr,"%lf",&x[i]);
    fclose(fptr);
    fptr = fopen("nodhn_sfc.out.orig", "w");
    for(i = 0; i < nod2D; i++) fprintf(fptr, "%lf\n", x[perm[i]]);
    fclose(fptr);
  }

  fptr = fopen("nodhn.out.tcarta.srtm.aceh_c", "r");
  if(fptr != NULL){
    for(i = 0; i < nod2D; i++) fscanf(fptr,"%lf",&x[i]);
    fclose(fptr);
    fptr = fopen("nodhn_sfc.out.tcarta.srtm.aceh_c", "w");
    for(i = 0; i < nod2D; i++) fprintf(fptr, "%lf\n", x[perm[i]]);
    fclose(fptr);
  }


/*   strcpy(&fpath[len], "/wef1"); */
/*   strcat(fpath, suffix); */
/*   fptr = fopen(fpath, "r"); */
/*   if(fptr != NULL){ */
/*     for(i = 0; i < nod2D; i++) */
/*       fscanf(fptr,"%lf",&x[i]); */
/*     fclose(fptr); */
/*     strcpy(&fpath[len], "/wef1_sfc"); */
/*     strcat(fpath, suffix); */
/*     fptr = fopen(fpath, "w"); */
/*     for(i = 0; i < nod2D; i++) */
/*       fprintf(fptr, "%lf\n", x[perm[i]]); */
/*     fclose(fptr); */
/*   } */
  free(x);
  free(perm);
  free(ind);


  
  strcpy(&fpath[len], "/elem2d");
  strcat(fpath, suffix);
  fptr = fopen(fpath, "r");
  if(fptr == NULL){
    printf("ERROR: %s konnte nicht geoeffnet werden.\n",fpath);
    return(1);
  }
  
  fscanf(fptr, "%d", &elem2D);
  el   = malloc(3*elem2D*sizeof(int));

  for(i = 0; i < elem2D; i++){
    fscanf(fptr, "%d %d %d", &j, &k, &l);
    /* Renumber the nodes with the new sorting */
    a = iperm[j-1];
    b = iperm[k-1];
    c = iperm[l-1];
    /* Sort each row numerically: smallest node number first */
    /* Do not change the orientation, though! */
    if (a < b && a < c) {
      el[3*i]   = a;
      el[3*i+1] = b;
      el[3*i+2] = c; }
    else if ( b < a &&  b < c){
      el[3*i]   = b;
      el[3*i+1] = c;
      el[3*i+2] = a;    }
    else {
      el[3*i]   = c;
      el[3*i+1] = a;
      el[3*i+2] = b;    }
  }
  fclose(fptr);

  /* Now, sort the rows numercally */
   qsort(el, elem2D, 3*sizeof(int), cmp_integer);

  
  strcpy(&fpath[len],"/elem2d_sfc");
  strcat(fpath, suffix);
  
  fptr = fopen(fpath, "w");
  fprintf(fptr, "%d\n", elem2D);
  for( i = 0; i < elem2D; i++) {
    fprintf(fptr, "%d %d %d\n", el[3*i], el[3*i+1] ,  el[3*i+2]  );
      }

  fclose(fptr);
 
  free(el);
  free(iperm);



}
 
 
  
int main(int argc, char **argv){
 
  clock_t t1, t2;
  if(argc < 2) {
    printf("Aufruf: %s path [suffix] (default: '.out')\n", argv[0]);
    return 1;
  }
  t1 = clock();
  if(argc < 3) 
    resort_grid(argv[1], ".out");
  else
    resort_grid(argv[1], argv[2]);
  
  t2 = clock();

  printf("benoetigte Zeit: %lf\n", (float) (t2-t1)/(float)CLOCKS_PER_SEC);
  return 0;
}
