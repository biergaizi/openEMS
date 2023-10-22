/*
 * Copyright (C) 2023 Yifeng Li (tomli@tomli.me)
 * Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ENGINE_KOKKOS_H
#define ENGINE_KOKKOS_H

#include "engine.h"
#include "operator_kokkos.h"

class Engine_Kokkos : public Engine
{
public:
	static Engine_Kokkos* New(const Operator_Kokkos* op);
	virtual ~Engine_Kokkos();

	virtual void Init();
	virtual void Reset();

	virtual unsigned int GetNumberOfTimesteps() {return numTS;};

	//this access functions muss be overloaded by any new engine using a different storage model
	inline virtual FDTD_FLOAT GetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		return (*volt_ptr)(n, x, y, z);
	}

	inline virtual FDTD_FLOAT GetVolt(unsigned int n, const unsigned int pos[3]) const
	{ 
		return (*volt_ptr)(n, pos[0], pos[1], pos[2]);
	}

	inline virtual FDTD_FLOAT GetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		return (*curr_ptr)(n, x, y, z);
	}

	inline virtual FDTD_FLOAT GetCurr(unsigned int n, const unsigned int pos[3]) const
	{
		return (*curr_ptr)(n, pos[0], pos[1], pos[2]);
	}

	inline virtual void SetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		(*volt_ptr)(n, x, y, z) = value;
	}

	inline virtual void SetVolt(unsigned int n, const unsigned int pos[3], FDTD_FLOAT value)
	{
		(*volt_ptr)(n, pos[0], pos[1], pos[2]) = value;
	}

	inline virtual void SetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		(*curr_ptr)(n, x, y, z) = value;
	}

	inline virtual void SetCurr(unsigned int n, const unsigned int pos[3], FDTD_FLOAT value)
	{
		(*curr_ptr)(n, pos[0], pos[1], pos[2]) = value;
	}

protected:
	Engine_Kokkos(const Operator_Kokkos* op);
	const Operator_Kokkos* Op;

	virtual void UpdateVoltages(unsigned int startX, unsigned int numX);
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);

public: //public access to the kokkos arrays for efficient extensions access... use careful...
	KokkosGlobalArray<float>* volt_ptr;
	KokkosGlobalArray<float>* curr_ptr;
};

#endif // ENGINE_KOKKOS_H
