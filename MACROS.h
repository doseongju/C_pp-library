#pragma once

#define ABS(x)      ((x) < 0 ? (-(x)) : (x))
#define SQR(x)      ((x) * (x)             )
#define CUBE(x)     ((x) * (x) * (x)       )
#define MAX(a,b)    ((a) > (b) ? (a) : (b) )
#define MAX3(a,b,c)	(MAX(MAX(a,b),c)       )
#define MAX4(a,b,c,d) (MAX(MAX(a,b), MAX(c,d)))
#define MIN(a,b)    ((a) < (b) ? (a) : (b) )
#define MIN3(a,b,c) (MIN(MIN(a,b),c)       )
#define MIN4(a,b,c,d) (MIN(MIN(a,b), MIN(c,d)))

#define POSPART(x) ((x) > 0 ? (x) : (0))
#define NEGPART(x) ((x) < 0 ? (x) : (0))

#define SIGN(x) ((x) > 0 ? 1 : (-1))

// sort a, b, c -> a1 < a2 < a3	
#define INCREASING_SORT3(a, b, c, a1, a2, a3)		if(a <= b){										\
														if(b <= c){a1 = a;a2 = b;a3 = c;}			\
														else if(a <= c){a1 = a;a2 = c;a3 = b;}		\
														else{a1 = c;a2 = a;a3 = b;}}				\
													 else{											\
													 if(a <= c){a1 = b;a2 = a;a3 = c;}				\
														else if(b <= c){a1 = b;a2 = c;a3 = a;}		\
														else{a1 = c;a2 = b;a3 = a;}}

#define ITERATION_GRID_1D(field, width) for(int i = field.i_start - width; i <= field.i_end + width; i++)
#define ITERATION_GRID_2D(field, width) for(int j = field.j_start - width; j <= field.j_end + width; j++)\
											for(int i = field.i_start - width; i <= field.i_end + width; i++)
#define ITERATION_GRID_2D_IJ(field, width) for(int i = field.i_start - width; i <= field.i_end + width; i++)\
											for(int j = field.j_start - width; j <= field.j_end + width; j++)
#define ITERATION_GRID_3D(field, width) for(int k = field.k_start - width; k <= field.k_end + width; k ++)\
										    for(int j = field.j_start - width; j <= field.j_end + width; j ++)\
												for(int i = field.i_start - width; i <= field.i_end + width; i ++)
#define LOOP_2D(i_start, i_end, j_start, j_end) for(int j = j_start; j <= j_end; j ++)\
													for(int i = i_start; i <= i_end; i ++)
#define PI (4. * atan(1.))

#define CLAMP(index, start, end) (((index) < (start)) ? (start) : ((index) > (end)) ? (end) : (index))

inline double MINMOD(const double& u, const double& v)
{
	if(u*v > 0)
		return SIGN(u) * MIN(ABS(u), ABS(v));
	else
		return 0;

}
#define MINMOD3(a, b, c) MINMOD(a, MINMOD(b, c))

#define SAFE_DELETE(pointer) if(pointer != 0){delete pointer; pointer=0;}
#define SAFE_DELETE_ARRAY(pointer) if(pointer != 0){delete [] pointer; pointer=0;}