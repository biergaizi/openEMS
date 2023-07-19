/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_EXT_EXCITATION_H
#define ENGINE_EXT_EXCITATION_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sycl.h"
#include "FDTD/operator.h"
#include <sycl/sycl.hpp>

class Operator_Ext_Excitation;

class Engine_Ext_Excitation : public Engine_Extension
{
public:
	Engine_Ext_Excitation(Operator_Ext_Excitation* op_ext);
	virtual ~Engine_Ext_Excitation();

	virtual void Apply2Voltages();
	virtual void Apply2Current();

	virtual void Apply2Voltages(int timestep, unsigned int start[3], unsigned int stop[3]);
	virtual void Apply2Current(int timestep, unsigned int start[3], unsigned int stop[3]);

	virtual void Apply2Voltages(sycl::queue Q);
	virtual void Apply2Current(sycl::queue Q);

	virtual void InitializeSYCL(sycl::queue Q);
	static void Apply2VoltagesSYCLKernel(
		Engine_sycl* eng,
		int n,
		int p,
		int numTS, int length,
		FDTD_FLOAT* exc_volt, int exc_pos,
		unsigned int* Volt_index[3],
		unsigned short* Volt_dir,
		FDTD_FLOAT* Volt_amp,
		unsigned int* Volt_delay
	);
	static void Apply2CurrentSYCLKernel(
		Engine_sycl* eng,
		int n,
		int p,
		int numTS, int length,
		FDTD_FLOAT* exc_curr, int exc_pos,
		unsigned int* Curr_index[3],
		unsigned short* Curr_dir,
		FDTD_FLOAT* Curr_amp,
		unsigned int* Curr_delay
	);

protected:
	Operator_Ext_Excitation* m_Op_Exc;

	bool InsideTile(
		unsigned int start[3], unsigned int stop[3],
		unsigned int ext_x, unsigned int ext_y, unsigned int ext_z
	);
};

#endif // ENGINE_EXT_EXCITATION_H
