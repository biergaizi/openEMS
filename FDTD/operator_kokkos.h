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

#ifndef OPERATOR_KOKKOS_H
#define OPERATOR_KOKKOS_H

#include "operator.h"
#include "tools/fdtd_grid3.hpp"

class Operator_Kokkos : public Operator
{
	friend class Engine_Interface_KOKKOS_FDTD;
public:
	//! Create a new operator
	static Operator_Kokkos* New();
	virtual ~Operator_Kokkos();

	virtual Engine* CreateEngine();

	inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		return (*vv_ptr)(n, x, y, z);
	}

	inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		return (*vi_ptr)(n, x, y, z);
	}

	inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		return (*ii_ptr)(n, x, y, z);
	}

	inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		return (*iv_ptr)(n, x, y, z);
	}

	inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		(*vv_ptr)(n, x, y, z) = value;
	}
	inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		(*vi_ptr)(n, x, y, z) = value;
	}
	inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		(*ii_ptr)(n, x, y, z) = value;
	}
	inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		(*iv_ptr)(n, x, y, z) = value;
	}

protected:
	//! use New() for creating a new Operator
	Operator_Kokkos();

	virtual void Init();
	void Delete();
	virtual void Reset();
	virtual void InitOperator();

	// engine/post-proc needs access
public:
	KokkosGrid* grid_ptr;
	KokkosGlobalArray<float>* vv_ptr;
	KokkosGlobalArray<float>* vi_ptr;
	KokkosGlobalArray<float>* ii_ptr;
	KokkosGlobalArray<float>* iv_ptr;
};

#endif // OPERATOR_KOKKOS_H
