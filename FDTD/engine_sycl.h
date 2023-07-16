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

#ifndef ENGINE_SYCL_H
#define ENGINE_SYCL_H

#include <sycl/sycl.hpp>
#include "engine.h"
#include "operator_sycl.h"

class Engine_sycl : public Engine
{
public:
	static Engine_sycl* New(const Operator_sycl* op);
	virtual ~Engine_sycl();

	virtual void Init();
	virtual void Reset();

	virtual unsigned int GetNumberOfTimesteps() {return numTS;};

	//this access functions muss be overloaded by any new engine using a different storage model
	inline virtual FDTD_FLOAT GetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const
	{
		SYCL_N_3DArray<sycl::float4> &f4_volt = *f4_volt_ptr;
		return f4_volt(n, x, y, z%numVectors).f[z/numVectors];
	}

	inline virtual FDTD_FLOAT GetVolt( unsigned int n, const unsigned int pos[3] ) const
	{
		SYCL_N_3DArray<sycl::float4> &f4_volt = *f4_volt_ptr;
		return f4_volt(n, pos[0], pos[1], pos[2]%numVectors).f[pos[2]/numVectors];
	}

	inline virtual FDTD_FLOAT GetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const
	{
		SYCL_N_3DArray<sycl::float4> &f4_curr = *f4_curr_ptr;
		return f4_curr(n, x, y, z%numVectors).f[z/numVectors];
	}

	inline virtual FDTD_FLOAT GetCurr( unsigned int n, const unsigned int pos[3] ) const
	{
		SYCL_N_3DArray<sycl::float4> &f4_curr = *f4_curr_ptr;
		return f4_curr(n, pos[0], pos[1], pos[2]%numVectors).f[pos[2]/numVectors];
	}

	inline virtual void SetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		SYCL_N_3DArray<sycl::float4> &f4_volt = *f4_volt_ptr;
		f4_volt(n, x, y, z%numVectors).f[z/numVectors]=value;
	}

	inline virtual void SetVolt( unsigned int n, const unsigned int pos[3], FDTD_FLOAT value )
	{
		SYCL_N_3DArray<sycl::float4> &f4_volt = *f4_volt_ptr;
		f4_volt(n, pos[0], pos[1], pos[2]%numVectors).f[pos[2]/numVectors]=value;
	}

	inline virtual void SetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		SYCL_N_3DArray<sycl::float4> &f4_curr = *f4_curr_ptr;
		f4_curr(n, x, y, z%numVectors).f[z/numVectors]=value;
	}

	inline virtual void SetCurr( unsigned int n, const unsigned int pos[3], FDTD_FLOAT value )
	{
		SYCL_N_3DArray<sycl::float4> &f4_curr = *f4_curr_ptr;
		f4_curr(n, pos[0], pos[1], pos[2]%numVectors).f[pos[2]/numVectors]=value;
	}

protected:
	Engine_sycl(const Operator_sycl* op);
	const Operator_sycl* Op;

	virtual void UpdateVoltages(unsigned int start[3], unsigned int stop[3]);
	virtual void UpdateCurrents(unsigned int start[3], unsigned int stop[3]);

	unsigned int numVectors;
	sycl::queue m_sycl_queue({sycl::property::queue::in_order()});

private:
	void Engine_sycl::UpdateVoltagesKernel(
	        const SYCL_N_3DArray<sycl::float4>& volt,
	        const SYCL_N_3DArray<sycl::float4>& curr,
	        const SYCL_N_3DArray<sycl::float4>& vv,
	        const SYCL_N_3DArray<sycl::float4>& vi,
	        int x, int y, int z
	);

	void Engine_sycl::UpdateCurrentsKernel(
	        const SYCL_N_3DArray<sycl::float4>& curr,
	        const SYCL_N_3DArray<sycl::float4>& volt,
	        const SYCL_N_3DArray<sycl::float4>& iv,
	        const SYCL_N_3DArray<sycl::float4>& ii,
	        int x, int y, int z
	);

public: //public access to the sycl arrays for efficient extensions access... use careful...
	SYCL_N_3DArray<sycl::float4>* f4_volt_ptr;
	SYCL_N_3DArray<sycl::float4>* f4_curr_ptr;
};

#endif // ENGINE_SYCL_H
