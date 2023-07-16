/*
*	Copyright (C) 2023 Yifeng Li (tomli@tomli.me)
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

#include "engine_sycl.h"

//! \brief construct an Engine_sycl instance
//! it's the responsibility of the caller to free the returned pointer
Engine_sycl* Engine_sycl::New(const Operator_sycl* op)
{
	cout << "Create FDTD engine (SYCL)" << endl;
	Engine_sycl* e = new Engine_sycl(op);
	e->Init();
	return e;
}

Engine_sycl::Engine_sycl(const Operator_sycl* op) : Engine(op)
{
	m_type = SYCL;
	Op = op;
	f4_volt_ptr = 0;
	f4_curr_ptr = 0;

	numVectors =  ceil((double)numLines[2]/4.0);
}

Engine_sycl::~Engine_sycl()
{
	Reset();
}

void Engine_sycl::Init()
{
	Engine::Init();

	Delete_N_3DArray(volt,numLines);
	volt=NULL; // not used
	Delete_N_3DArray(curr,numLines);
	curr=NULL; // not used

	f4_volt_ptr = Create_SYCL_N_3DArray<f4vector>(numLines);
	f4_curr_ptr = Create_SYCL_N_3DArray<f4vector>(numLines);
}

void Engine_sycl::Reset()
{
	Engine::Reset();
	Delete_SYCL_N_3DArray(f4_volt_ptr,numLines);
	f4_volt_ptr = 0;
	Delete_SYCL_N_3DArray(f4_curr_ptr,numLines);
	f4_curr_ptr = 0;
}

/*
 * UpdateVoltages
 *
 * Submit the UpdateVoltagesKernel to the GPU for execution.
 */
void Engine_sycl::UpdateVoltages(
	unsigned int start[3],
	unsigned int stop[3]
)
{
	m_sycl_queue.submit([&](sycl::handler &h)
	{
		h.parallel_for<class Voltage>(
			sycl::range(
				stop[0] - start[0],
				stop[1] - start[1],
				stop[2] - start[2]
			),
			[=](sycl::item<3> itm)
			{
				/* this C++ lambda is the function body for GPU */
				int x = itm.get_id(0) + start[0];
				int y = itm.get_id(1) + start[1];
				int z = itm.get_id(2) + start[2];

				UpdateVoltages(volt, curr, vv, vi, x, y, z);
			}
		);
	});

	/* block until the GPU finishes */
	m_sycl_queue.wait();
}

/*
 * UpdateVoltagesKernel
 *
 * Calculate new electric field array "volt" based on
 * magnetic field "curr" and two electromagnetic field
 * operators "vv" and "vi", precalculated before starting
 * up simulation.
 *
 * Note: This is a data-parallel SYCL kernel. Many copies of
 * UpdateVoltages() are simutaneously executed on the GPU at
 * different cells in the 3D space with different "tid".
 * This is called Single Program Multiple Data (SPMD).
 */
void Engine_sycl::UpdateVoltagesKernel(
        const SYCL_N_3DArray<sycl::float4>& volt,
        const SYCL_N_3DArray<sycl::float4>& curr,
        const SYCL_N_3DArray<sycl::float4>& vv,
        const SYCL_N_3DArray<sycl::float4>& vi,
        int x, int y, int z
)
{
        // note: each (x, y, z) cell has three polarizations
        // x, y, z, these are different from the cell's
        // coordinates (x, y, z)

        //for x polarization
        sycl::float4 volt0 = volt(0, x, y, z);
        volt0 *= vv(0, x, y, z);
        volt0 +=
            vi(0, x, y, z) * (
                curr(2, x, y       , z       ) -
                curr(2, x, y-prev_y, z       ) -
                curr(1, x, y       , z       ) +
                curr(1, x, y       , z-prev_z)
            );

        //for y polarization
        sycl::float4 volt1 = volt(1, x, y, z);
        volt1 *= vv(1, x, y, z);
        volt1 +=
            vi(1, x, y, z) * (
                curr(0, x       , y, z       ) -
                curr(0, x       , y, z-prev_z) -
                curr(2, x       , y, z       ) +
                curr(2, x-prev_x, y, z       )
            );

        //for z polarization
        sycl::float4 volt2 = volt(2, x, y, z);
        volt2 *= vv(2, x, y, z);
        volt2 +=
            vi(2, x, y, z) * (
                curr(1, x       , y       , z) -
                curr(1, x-prev_x, y       , z) -
                curr(0, x       , y       , z) +
                curr(0, x       , y-prev_y, z)
            );

        volt(0, x, y, z) = volt0;
        volt(1, x, y, z) = volt1;
        volt(2, x, y, z) = volt2;
}

/*
 * UpdateCurrents
 *
 * Submit the UpdateCurrentsKernel to the GPU for execution.
 */
void Engine_sycl::UpdateCurrents(
	unsigned int start[3],
	unsigned int stop[3]
)
{
	m_sycl_queue.submit([&](sycl::handler &h)
	{
		h.parallel_for<class Voltage>(
			sycl::range(
				stop[0] - start[0],
				stop[1] - start[1],
				stop[2] - start[2]
			),
			[=](sycl::item<3> itm)
			{
				/* this C++ lambda is the function body for GPU */
				int x = itm.get_id(0) + start[0];
				int y = itm.get_id(1) + start[1];
				int z = itm.get_id(2) + start[2];

				UpdateCurrents(volt, curr, vv, vi, x, y, z);
			}
		);
	});

	/* block until the GPU finishes */
	m_sycl_queue.wait();
}

/*
 * UpdateCurrentsKernel
 *
 * Calculate new magnetic field array "curr" based on
 * electric field "volt" and two electromagnetic field
 * operators "ii" and "iv", precalculated before starting
 * simulation.
 *
 * Note: This is a data-parallel SYCL kernel. Many copies of
 * UpdateCurrents() are simutaneously executed on the GPU at
 * different cells in the 3D space with different "tid".
 * This is called Single Program Multiple Data (SPMD).
 */
void Engine_sycl::UpdateCurrentsKernel(
        const SYCL_N_3DArray<sycl::float4>& curr,
        const SYCL_N_3DArray<sycl::float4>& volt,
        const SYCL_N_3DArray<sycl::float4>& iv,
        const SYCL_N_3DArray<sycl::float4>& ii,
        int x, int y, int z
)
{
        // note: each (x, y, z) cell has three polarizations
        // x, y, z, these are different from the cell's
        // coordinates (x, y, z)

        //for x polarization
        sycl::float4 curr0 = curr(0, x, y, z);
        curr0 *= ii(0, x, y, z);
        curr0 +=
            iv(0, x, y, z) * (
                volt(2, x, y  , z  ) -
                volt(2, x, y+1, z  ) -
                volt(1, x, y  , z  ) +
                volt(1, x, y  , z+1)
            );

        //for y polarization
        sycl::float4 curr1 = curr(1, x, y, z);
        curr1 *= ii(1, x, y, z);
        curr1 +=
            iv(1, x, y, z) * (
                volt(0, x  , y, z  ) -
                volt(0, x  , y, z+1) -
                volt(2, x  , y, z  ) +
                volt(2, x+1, y, z  )
            );

        //for z polarization
        sycl::float4 curr2 = curr(2, x, y, z);
        curr2 *= ii(2, x, y, z);
        curr2 +=
            iv(2, x, y, z) * (
                volt(1, x  , y  , z) -
                volt(1, x+1, y  , z) -
                volt(0, x  , y  , z) +
                volt(0, x  , y+1, z)
            );

        curr(0, x, y, z) = curr0;
        curr(1, x, y, z) = curr1;
        curr(2, x, y, z) = curr2;
}
