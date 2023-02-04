/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "array_ops.h"
#include <ostream>
#include <string.h>

using namespace std;

#ifdef WIN32
#include <malloc.h>
#define MEMALIGN( array, alignment, size ) !(*array = _mm_malloc( size, alignment ))
#define FREE( array ) _mm_free( array )
#else
#define MEMALIGN( array, alignment, size ) posix_memalign( array, alignment, size )
#define FREE( array ) free( array )
#endif

void Delete1DArray_v4sf(f4vector* array)
{
	if (array==NULL) return;
	FREE( array );
}


void Delete3DArray_v4sf(f4vector*** array, const unsigned int* numLines)
{
	if (array==NULL) return;
	unsigned int pos[3];
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			FREE( array[pos[0]][pos[1]] );
			//delete[] array[pos[0]][pos[1]];
		}
		FREE( array[pos[0]] );
		//delete[] array[pos[0]];
	}
	FREE( array );
	//delete[] array;
}

void Delete_N_3DArray_v4sf(f4vector**** array, const unsigned int* numLines)
{
	if (array==NULL) return;
	for (int n=0; n<3; ++n)
	{
		Delete3DArray_v4sf(array[n],numLines);
	}
	FREE( array );
	//delete[] array;
}

f4vector* Create1DArray_v4sf(const unsigned int numLines)
{
	f4vector* array=NULL;
	if (MEMALIGN( (void**)&array, 16, F4VECTOR_SIZE*numLines ))
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	for (unsigned int pos=0; pos<numLines; ++pos)
	{
		array[pos].f[0] = 0;
		array[pos].f[1] = 0;
		array[pos].f[2] = 0;
		array[pos].f[3] = 0;
	}
	return array;
}

//! \brief this function allocates a 3D array, which is aligned to 16 byte
f4vector*** Create3DArray_v4sf(const unsigned int* numLines)
{
	unsigned int numZ = ceil((double)numLines[2]/4.0);

	f4vector*** array=NULL;
	unsigned int pos[3];
	if (MEMALIGN( (void**)&array, 16, F4VECTOR_SIZE*numLines[0] ))
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	//array = new f4vector**[numLines[0]];
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		if (MEMALIGN( (void**)&array[pos[0]], 16, F4VECTOR_SIZE*numLines[1] ))
		{
			cerr << "cannot allocate aligned memory" << endl;
			exit(3);
		}
		//array[pos[0]] = new f4vector*[numLines[1]];
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			if (MEMALIGN( (void**)&array[pos[0]][pos[1]], 16, F4VECTOR_SIZE*numZ ))
			{
				cerr << "cannot allocate aligned memory" << endl;
				exit(3);
			}
			//array[pos[0]][pos[1]] = new f4vector[numZ];
			for (pos[2]=0; pos[2]<numZ; ++pos[2])
			{
				array[pos[0]][pos[1]][pos[2]].f[0] = 0;
				array[pos[0]][pos[1]][pos[2]].f[1] = 0;
				array[pos[0]][pos[1]][pos[2]].f[2] = 0;
				array[pos[0]][pos[1]][pos[2]].f[3] = 0;
			}
		}
	}
	return array;
}

N_3DArray_v4sf *Create_N_3DArray_Flat_v4sf(const unsigned int* numLines)
{
	N_3DArray_v4sf *n_3d_array_v4sf;
	unsigned int n_max = 3;
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];
	unsigned int z_max = ceil((double)numLines[2] / 4.0);

	// Size of the header itself.
	// If the definition of N_3DArray_v4sf has been changed,
	// all its data type must be a multiple of 16 (F4VECTOR_SIZE).
	static_assert(sizeof(N_3DArray_v4sf) % F4VECTOR_SIZE == 0);
	size_t size = sizeof(N_3DArray_v4sf);

	// and the actual memory of the array[1] flexible array member
	size += F4VECTOR_SIZE * n_max * x_max * y_max * z_max;

	// array[0] is counted twice, so remove one element.
	size -= F4VECTOR_SIZE;

	if (MEMALIGN( (void**)&n_3d_array_v4sf, 16, size))
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	memset(n_3d_array_v4sf, 0, size);

	//n_3d_array_v4sf->n_max = n_max;
	//n_3d_array_v4sf->x_max = x_max;
	//n_3d_array_v4sf->y_max = y_max;
	//n_3d_array_v4sf->z_max = z_max;
	n_3d_array_v4sf->x_stride = y_max * z_max * n_max;
	n_3d_array_v4sf->y_stride = z_max * n_max;

	return n_3d_array_v4sf;
}

N_3DArray* Create_N_3DArray_Flat(const unsigned int* numLines)
{
	N_3DArray *n_3d_array;
	unsigned int n_max = 3;
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];
	unsigned int z_max = numLines[2];

	size_t size = sizeof(N_3DArray);

	// and the actual memory of the array[1] flexible array member
	size += sizeof(float) * n_max * x_max * y_max * z_max;

	// array[0] is counted twice, so remove one element.
	size -= sizeof(float);

	if (MEMALIGN( (void**)&n_3d_array, 16, size))
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	memset(n_3d_array, 0, size);

	n_3d_array->x_stride = y_max * z_max * n_max;
	n_3d_array->y_stride = z_max * n_max;

	return n_3d_array;
}

void Delete_N_3DArray_Flat(N_3DArray* array, const unsigned int* numLines)
{
	free(array);
}
