#ifndef FLAT_ARRAY_OPS_H
#define FLAT_ARRAY_OPS_H

#include <cstring>
#include "array_ops.h"
#include <ostream>

//***************************
// Flat2DArray
//***************************

template <typename T>
struct Flat2DArray
{
	inline T& operator() (const unsigned int x, const unsigned int y)
	{
		return array[x * x_stride + y];
	}

	inline T operator() (const unsigned int &x, const unsigned int &y) const
	{
		return array[x * x_stride + y];
	}

	unsigned long x_stride, pad;

	// This is a flexible array member, the actual size would be
	// determined by the actual size used to call malloc (actually
	// posix_memalign() if the type is a SIMD f4vector).
	//
	// It's technically a standard-nonconforming undefined behavior,
	// but is a well-known technique and it's important here to
	// avoid the cost of extra dereferencing.
	T array[1];
};

template <typename T>
Flat2DArray<T>* CreateFlat2DArray(const unsigned int* numLines)
{
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];

	// size of the Flat_2DArray class itself.
	size_t size = sizeof(Flat2DArray<T>);

	// and the actual memory of the array[1] flexible array member
	size += sizeof(T) * x_max * y_max;

	// array[0] is counted twice, so remove one element.
	size -= sizeof(T);

	void *buf;
	if (MEMALIGN(&buf, 16, size))
	{
		std::cerr << "cannot allocate aligned memory" << std::endl;
		exit(3);
	}
	memset(buf, 0, size);

	// create Flat_2DArray inside manually allocated memory "buf".
	Flat2DArray<T>* array = new (buf) Flat2DArray<T>;
	array->x_stride = y_max;

	return array;
}

template <typename T>
void DeleteFlat2DArray(Flat2DArray<T>* array, const unsigned int* numLines)
{
	free(array);
}

//***************************
// Flat_N_3DArray
//***************************

template <typename T>
struct Flat_N_3DArray
{
	inline T& operator() (const unsigned int n, const unsigned int x, const unsigned int y, const unsigned int z)
	{
		return array[
		           x * x_stride +
		           y * y_stride +
		           z * (3) +
		           n
		       ];
	}

	inline T operator() (const unsigned int &n, const unsigned int &x, const unsigned int &y, const unsigned int &z) const
	{
		return array[
		           x * x_stride +
		           y * y_stride +
		           z * (3) +
		           n
		       ];
	}

	unsigned long x_stride, y_stride;

	// This is a flexible array member, the actual size would be
	// determined by the actual size used to call malloc (actually
	// posix_memalign() if the type is a SIMD f4vector).
	//
	// It's technically a standard-nonconforming undefined behavior,
	// but is a well-known technique and it's important here to
	// avoid the cost of extra dereferencing.
	T array[1];
};

template <typename T>
Flat_N_3DArray<T>* Create_Flat_N_3DArray(const unsigned int* numLines)
{
	unsigned int n_max = 3;
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];
	unsigned int z_max = numLines[2];

	// size of the Flat_N_3DArray class itself.
	size_t size = sizeof(Flat_N_3DArray<T>);

	// and the actual memory of the array[1] flexible array member
	size += sizeof(T) * n_max * x_max * y_max * z_max;

	// array[0] is counted twice, so remove one element.
	size -= sizeof(T);

	void *buf;
	if (MEMALIGN(&buf, 16, size))
	{
		std::cerr << "cannot allocate aligned memory" << std::endl;
		exit(3);
	}
	memset(buf, 0, size);

	// create Flat_N_3DArray inside manually allocated memory "buf".
	Flat_N_3DArray<T>* array = new (buf) Flat_N_3DArray<T>;
	array->x_stride = y_max * z_max * n_max;
	array->y_stride = z_max * n_max;

	return array;
}

template <>
inline Flat_N_3DArray<f4vector>* Create_Flat_N_3DArray<f4vector>(const unsigned int* numLines)
{
	unsigned int n_max = 3;
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];
	unsigned int z_max = ceil((double) numLines[2] / 4.0);

	// Size of the header itself.
	// If the definition of Flat_N_3DArray<f4vector> has been changed,
	// all its data type must be a multiple of 16 (F4VECTOR_SIZE).
	static_assert(sizeof(Flat_N_3DArray<f4vector>) % F4VECTOR_SIZE == 0);
	size_t size = sizeof(Flat_N_3DArray<f4vector>);

	// and the actual memory of the array[1] flexible array member
	size += F4VECTOR_SIZE * n_max * x_max * y_max * z_max;

	// array[0] is counted twice, so remove one element.
	size -= F4VECTOR_SIZE;

	printf("allocating %llu bytes\n", size);

	void *buf;
	if (MEMALIGN(&buf, 32, size))
	{
		std::cerr << "cannot allocate aligned memory" << std::endl;
		exit(3);
	}
	memset(buf, 0, size);

	// create Flat_N_3DArray inside manually allocated memory "buf".
	Flat_N_3DArray<f4vector>* array = new (buf) Flat_N_3DArray<f4vector>;
	array->x_stride = y_max * z_max * n_max;
	array->y_stride = z_max * n_max;
	printf("obj: %p, x_stride: %p, y_stride: %p, array: %p\n",
	       array, &array->x_stride, &array->y_stride, array->array);

	return array;
}

template <typename T>
void Delete_Flat_N_3DArray(Flat_N_3DArray<T>* array, const unsigned int* numLines)
{
	free(array);
}

template <>
inline void Delete_Flat_N_3DArray<f4vector>(Flat_N_3DArray<f4vector>* array, const unsigned int* numLines)
{
	FREE(array);
}

#endif
