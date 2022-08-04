
// Stephan Frickenhaus, Alfred Wegener Institute (2022), based on 

//Rakowsky N. & Fuchs A. Efficient local resorting techniques with space filling curves applied
//to the tsunami simulation model TsunAWI. In IMUM 2011 - The 10th International Workshop
//on Multiscale (Un-)structured Mesh Numerical Modelling for coastal, shelf and global ocean
//  dynamics, August 2011. http://hdl.handle.net/10013/epic.39576.d001

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include<Rcpp.h>

enum sfc {A, D, C, U};

int cmp_integer(const void *a, const void *b){
  int *i,*j;
  i = (int *)a;
  j = (int *)b;
  return (i[0]-j[0]);
}

using namespace Rcpp;
//[[Rcpp::export]]
IntegerVector mysfc(NumericVector lat, NumericVector lon) 
{ 
int dim=lat.size();
IntegerVector p(dim,0);
if (lon.size()!=lat.size()) { stop("x and y differ in length!");  }

  int i, j;
  int maxlevel = 14;
  int *indx, indx_add;
//   int halve;
  int *perm;
  
  int lc;
  double xmin, xmax,  ymin, ymax, mx, ny, xhalve, yhalve;
  // double xscale, yscale;
  int nmax;
  enum sfc c;

  nmax = 2;
  for (i=0; i<maxlevel; i++) nmax *= 2;
 // Rprintf("dim %i\n",dim);
  
  xmin = xmax = lat[0];
  ymin = ymax = lon[0];
  
  for(i = 1; i < dim; i++){
    xmin = (lat[i] < xmin) ? lat[i] : xmin;
    xmax = (lat[i] > xmax) ? lat[i] : xmax;
    ymin = (lon[i] < ymin) ? lon[i] : ymin;
    ymax = (lon[i] > ymax) ? lon[i] : ymax;
  }

 // xscale = ((double) nmax)/(xmax-xmin);
  // yscale = ((double) nmax)/(ymax-ymin);
  
  indx = (int *) calloc(dim, sizeof(int));
  
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
	  REprintf("ERROR: falsches lc\n");	  
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
	  REprintf("ERROR: falsches lc\n");	  
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
	  REprintf("ERROR: falsches lc\n");	  
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
	  REprintf("ERROR: falsches lc\n");	  
	  break;
	  
	}	  
	break;	
	  
      }	
    }    
  }
   
  perm = (int *) malloc(2*dim*sizeof(int));
  for(i = 0; i < dim; i++){
    perm[2*i] = indx[i];
    perm[2*i+1] = i;
  }  
  qsort(perm, dim, 2*sizeof(int), cmp_integer);
  for(i = 0; i < dim; i++)
    p[i] = perm[2*i+1] + 1; // return R-index

  free(indx);
  free(perm);
//  
  return(p);
}



