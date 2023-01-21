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

#include "engine_interface_sse_fdtd.h"

#define x_max m_Eng_SSE->Engine_sse::x_max
#define y_max m_Eng_SSE->Engine_sse::y_max
#define z_max m_Eng_SSE->Engine_sse::z_max
#define f4_volt(n, x, y, z)			\
	(m_Eng_SSE->Engine_sse::_f4_volt	\
		[				\
		 n * (x_max * y_max * z_max) +  \
		 x * (y_max * z_max) +		\
		 y * (z_max) +			\
		 z				\
		]				\
	)

#define f4_curr(n, x, y, z)			\
	(m_Eng_SSE->Engine_sse::_f4_curr	\
		[				\
		 n * (x_max * y_max * z_max) +  \
		 x * (y_max * z_max) +		\
		 y * (z_max) +			\
		 z				\
		]				\
	)

Engine_Interface_SSE_FDTD::Engine_Interface_SSE_FDTD(Operator_sse* op) : Engine_Interface_FDTD(op)
{
	m_Op_SSE = op;
	m_Eng_SSE = dynamic_cast<Engine_sse*>(m_Op_SSE->GetEngine());
	if (m_Eng_SSE==NULL)
	{
		cerr << "Engine_Interface_SSE_FDTD::Engine_Interface_SSE_FDTD: Error: SSE-Engine is not set! Exit!" << endl;
		exit(1);
	}
}

Engine_Interface_SSE_FDTD::~Engine_Interface_SSE_FDTD()
{
	m_Op_SSE=NULL;
	m_Eng_SSE=NULL;
}

double Engine_Interface_SSE_FDTD::CalcFastEnergy() const
{
	f4vector E_energy;
	E_energy.f[0]=0;
	E_energy.f[1]=0;
	E_energy.f[2]=0;
	E_energy.f[3]=0;
	f4vector H_energy;
	H_energy = E_energy;

	if (m_Eng_SSE->GetType()!=Engine::SSE)
		return Engine_Interface_FDTD::CalcFastEnergy();

	unsigned int pos[3];
	for (pos[0]=0; pos[0]<m_Op_SSE->GetNumberOfLines(0)-1; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_Op_SSE->GetNumberOfLines(1)-1; ++pos[1])
		{
			for (pos[2]=0; pos[2]<m_Op_SSE->numVectors; ++pos[2])
			{
				E_energy.v += f4_volt(0, pos[0], pos[1], pos[2]).v * f4_volt(0, pos[0], pos[1], pos[2]).v;
				E_energy.v += f4_volt(1, pos[0], pos[1], pos[2]).v * f4_volt(1, pos[0], pos[1], pos[2]).v;
				E_energy.v += f4_volt(2, pos[0], pos[1], pos[2]).v * f4_volt(2, pos[0], pos[1], pos[2]).v;

				H_energy.v += f4_curr(0, pos[0], pos[1], pos[2]).v * f4_curr(0, pos[0], pos[1], pos[2]).v;
				H_energy.v += f4_curr(1, pos[0], pos[1], pos[2]).v * f4_curr(1, pos[0], pos[1], pos[2]).v;
				H_energy.v += f4_curr(2, pos[0], pos[1], pos[2]).v * f4_curr(2, pos[0], pos[1], pos[2]).v;
			}
		}
	}

	return __EPS0__*(E_energy.f[0]+E_energy.f[1]+E_energy.f[2]+E_energy.f[3]) + __MUE0__*(H_energy.f[0]+H_energy.f[1]+H_energy.f[2]+H_energy.f[3]);
}
