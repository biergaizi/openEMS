#ifndef FLAT_ARRAY_OPS_H
#define FLAT_ARRAY_OPS_H

#include <cstring>
#include "array_ops.h"
#include <ostream>

template <typename T>
struct Flat_N_3DArray
{
	inline T& operator() (const unsigned int n, const unsigned int x, const unsigned int y, const unsigned int z)
	{
		return array[
		           n * n_stride +
		           x * x_stride +
		           y * y_stride +
		           z
		       ];
	}

	inline const T operator() (const unsigned int &n, const unsigned int &x, const unsigned int &y, const unsigned int &z) const
	{
		return array[
		           n * n_stride +
		           x * x_stride +
		           y * y_stride +
		           z
		       ];
	}

	unsigned long n_stride, x_stride, y_stride;

	// This is a flexible array member, the size of 1 is just a place-
	// holder and its actual size is determined by how much memory is
	// allocated beyond the end of the array at runtime via malloc()
	// or posix_memalign().
	//
	// It's technically a standard-nonconforming undefined behavior,
	// but is a well-known technique and it's important here to
	// avoid the cost of extra dereferencing. GCC officially supports
	// the "size 0" extension, so we use this well-defined option when
	// it's available. C99 standardized "array[]" for the same purpose
	// but it's never added to C++.
#ifdef __GNUC__
	T array[0];
#else
	T array[1];
#endif
};

template <typename T>
Flat_N_3DArray<T>* Create_Flat_N_3DArray(const unsigned int* numLines)
{
	unsigned int n_max = 3;
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];
	unsigned int z_max = numLines[2];

	unsigned long n_stride = x_max * y_max * z_max;
	unsigned long x_stride = y_max * z_max;
	unsigned long y_stride = z_max;

	if (n_stride % 64 != 0) {
		n_stride += 64 - (n_stride % 64);
	}

	// size of the Flat_N_3DArray class itself.
	size_t size = sizeof(Flat_N_3DArray<T>);

	// and the actual memory of the array[0] flexible array member
	size += sizeof(T) * n_stride * n_max;

	// then allocate a chunk of memory with this size
	void *buf;
	if (MEMALIGN(&buf, 64, size))
	{
		std::cerr << "cannot allocate aligned memory" << std::endl;
		exit(3);
	}
	memset(buf, 0, size);

	// create Flat_N_3DArray inside manually allocated memory "buf".
	Flat_N_3DArray<T>* array = new (buf) Flat_N_3DArray<T>;
	array->n_stride = n_stride;
	array->x_stride = x_stride;
	array->y_stride = y_stride;

	return array;
}

template <>
inline Flat_N_3DArray<f4vector>* Create_Flat_N_3DArray<f4vector>(const unsigned int* numLines)
{
	unsigned int n_max = 3;
	unsigned int x_max = numLines[0];
	unsigned int y_max = numLines[1];
	unsigned int z_max = ceil((double) numLines[2] / 4.0);

	unsigned long n_stride = x_max * y_max * z_max;
	unsigned long x_stride = y_max * z_max;
	unsigned long y_stride = z_max;

	if (n_stride % 64 != 0) {
		n_stride += 64 - (n_stride % 64);
	}

	// If the definition of Flat_N_3DArray<f4vector> has been changed,
	// the underlying array must start at an address that is a multiple
	// of 16 (F4VECTOR_SIZE) to satisfy SSE alignment requirement. Thus,
	// adding new fields to the header may require manual padding.
	static_assert(offsetof(Flat_N_3DArray<f4vector>, array) % F4VECTOR_SIZE == 0);

	// Size of the header itself.
	size_t size = sizeof(Flat_N_3DArray<f4vector>);

	// and the actual memory of the array[0] flexible array member
	size += F4VECTOR_SIZE * n_stride * n_max;

	// then allocate a chunk of memory with this size
	void *buf;
	if (MEMALIGN(&buf, 64, size))
	{
		std::cerr << "cannot allocate aligned memory" << std::endl;
		exit(3);
	}
	memset(buf, 0, size);

	// create Flat_N_3DArray inside manually allocated memory "buf".
	Flat_N_3DArray<f4vector>* array = new (buf) Flat_N_3DArray<f4vector>;
	array->n_stride = n_stride;
	array->x_stride = x_stride;
	array->y_stride = y_stride;

	return array;
}

template <typename T>
void Delete_Flat_N_3DArray(Flat_N_3DArray<T>* array, const unsigned int* numLines)
{
	if (!array)
	{
		return;
	}

	FREE(array);
}

#endif
