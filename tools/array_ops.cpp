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
#if 0
	for (int n=0; n<3; ++n)
	{
		Delete3DArray_v4sf(array[n],numLines);
	}
	FREE( array );
	//delete[] array;
#endif
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

// Helper function to help calculating the offset of a
// contiguous iliffe vector, as if malloc is used - because
// the explicit formula is a brain-teaser.
//
// "*ptr" must be always pointed to the beginning of the
// underlying iliffe vector block, already returned by
// malloc().
static size_t ialloc_offset = 0;
static void *ialloc(void *ptr, size_t size)
{
	void *retval = (void *) (((char *) ptr) + ialloc_offset);
	ialloc_offset += size;
	return retval;
}

f4vector**** Create_N_3DArray_v4sf(const unsigned int* numLines)
{
	unsigned int n_max = 3;
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];
	unsigned int z_max = ceil((double)numLines[2] / 4.0);

	// create the N-3D array itself
	f4vector* array = NULL;
	size_t size = sizeof(f4vector) * n_max * x_max * y_max * z_max;
	if (MEMALIGN( (void**)&array, 16, size ) )
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	memset(array, 0, size);

	// create an iliffe vector to the N-3D array to allow
	// array[n][x][y][z] access
	void *iliffe_vector = malloc(sizeof(f4vector *) * (n_max +
							   n_max * x_max +
							   n_max * x_max * y_max));
	ialloc_offset = 0;

	f4vector**** array_ptr = NULL;
	array_ptr = (f4vector****) ialloc(iliffe_vector, sizeof(f4vector ***) * n_max);
	for (unsigned int n = 0; n < n_max; n++) {
		array_ptr[n] = (f4vector***) ialloc(iliffe_vector, sizeof(f4vector **) * x_max);
		for (unsigned int x = 0; x < x_max; x++) {
			array_ptr[n][x] = (f4vector **) ialloc(iliffe_vector, sizeof(f4vector *) * y_max);
			for (unsigned int y = 0; y < y_max; y++) {
				size_t offset = n * (x_max * y_max * z_max) +
						x * (y_max * z_max) +
						y * (z_max);
				array_ptr[n][x][y] = &array[offset];
			}
		}
	}
	return array_ptr;
}

