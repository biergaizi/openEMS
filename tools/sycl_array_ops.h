/*
*	Copyright (C) 2023 Yifeng Li <tomli@tomli.me>
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

#ifndef SYCL_ARRAY_OPS_H
#define SYCL_ARRAY_OPS_H

#define hipMemAdviseSetPreferredLocation 3
#define hipMemAdviseSetCoarseGrain 100

#include <cstring>
#include <ostream>

template <typename T>
struct SYCL_N_3DArray
{
        inline T& operator() (const unsigned int n, const unsigned int x, const unsigned int y, const unsigned int z) const
        {
                return array[
                           n * n_stride +
                           x * x_stride +
                           y * y_stride +
                           z
                       ];
        }

        unsigned long n_stride, x_stride, y_stride, size;
        T *array;
};

template <typename T>
SYCL_N_3DArray<T>* Create_SYCL_N_3DArray(sycl::queue Q, const unsigned int* numLines)
{
        unsigned int n_max = 3;
        unsigned int x_max = numLines[0];
        unsigned int y_max = numLines[1];
        unsigned int z_max = numLines[2];

        unsigned long n_stride = x_max * y_max * z_max;
        unsigned long x_stride = y_max * z_max;
        unsigned long y_stride = z_max;

        if (n_stride % 128 != 0)
        {
                n_stride += 128 - (n_stride % 128);
        }

        // allocate 1D linear buffer
        size_t size = n_stride * n_max;

        T *buf = sycl::malloc_shared<T>(size, Q);
	sycl::mem_advise(buf, size * sizeof(T), hipMemAdviseSetPreferredLocation, Q);
	sycl::mem_advise(buf, size * sizeof(T), hipMemAdviseSetCoarseGrain, Q);

	Q.submit([&](sycl::handler& h) {
		h.memset(buf, 0, size * sizeof(T));
	});
	Q.wait();

        // allocate wrapper class
        SYCL_N_3DArray<T>* array = new SYCL_N_3DArray<T>();
        array->n_stride = n_stride;
        array->x_stride = x_stride;
        array->y_stride = y_stride;
        array->size = size * sizeof(T);
        array->array = buf;

        return array;
}

template <>
inline SYCL_N_3DArray<sycl::float4>* Create_SYCL_N_3DArray<sycl::float4>(sycl::queue Q, const unsigned int* numLines)
{
        unsigned int n_max = 3;
        unsigned int x_max = numLines[0];
        unsigned int y_max = numLines[1];
        unsigned int z_max = ceil((double) numLines[2] / 4.0);

        unsigned long n_stride = x_max * y_max * z_max;
        unsigned long x_stride = y_max * z_max;
        unsigned long y_stride = z_max;

        if (n_stride % 128 != 0)
        {
                n_stride += 128 - (n_stride % 128);
        }

        // allocate 1D linear buffer
        size_t size = n_stride * n_max;

        sycl::float4* buf = sycl::malloc_shared<sycl::float4>(size, Q);
	sycl::mem_advise(buf, size * sizeof(sycl::float4), hipMemAdviseSetPreferredLocation, Q);
	sycl::mem_advise(buf, size * sizeof(sycl::float4), hipMemAdviseSetCoarseGrain, Q);

	Q.submit([&](sycl::handler& h) {
		h.memset(buf, 0, size * sizeof(sycl::float4));
	});
	Q.wait();

        // allocate wrapper class
        SYCL_N_3DArray<sycl::float4>* array = new SYCL_N_3DArray<sycl::float4>();
        array->n_stride = n_stride;
        array->x_stride = x_stride;
        array->y_stride = y_stride;
        array->size = size * sizeof(sycl::float4);
        array->array = buf;

        return array;
}

template <typename T>
void Delete_SYCL_N_3DArray(SYCL_N_3DArray<T>* array, const unsigned int* numLines)
{
	if (!array)
	{
		return;
	}
}

#endif
