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

#include "engine_ext_excitation.h"
#include "operator_ext_excitation.h"
#include "FDTD/engine_sse.h"

Engine_Ext_Excitation::Engine_Ext_Excitation(Operator_Ext_Excitation* op_ext) : Engine_Extension(op_ext)
{
	m_Op_Exc = op_ext;
	m_Priority = ENG_EXT_PRIO_EXCITATION;

	// this extension support the tiling engine
	m_TilingSupported = true;
}

Engine_Ext_Excitation::~Engine_Ext_Excitation()
{

}

void Engine_Ext_Excitation::Apply2Voltages()
{
	//soft voltage excitation here (E-field excite)
	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	int numTS = m_Eng->GetNumberOfTimesteps();
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_volt =  m_Op_Exc->m_Exc->GetVoltageSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	//switch for different engine types to access faster inline engine functions
	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int n=0; n<m_Op_Exc->Volt_Count; ++n)
			{
				exc_pos = numTS - (int)m_Op_Exc->Volt_delay[n];
				exc_pos *= (exc_pos>0);
				exc_pos %= p;
				exc_pos *= (exc_pos<(int)length);
				ny = m_Op_Exc->Volt_dir[n];
				pos[0]=m_Op_Exc->Volt_index[0][n];
				pos[1]=m_Op_Exc->Volt_index[1][n];
				pos[2]=m_Op_Exc->Volt_index[2][n];
				m_Eng->Engine::SetVolt(ny,pos, m_Eng->Engine::GetVolt(ny,pos) + m_Op_Exc->Volt_amp[n]*exc_volt[exc_pos]);
			}
			break;
		}
	case Engine::SSE:
		{
			for (unsigned int n=0; n<m_Op_Exc->Volt_Count; ++n)
			{
				Engine_sse* eng_sse = (Engine_sse*) m_Eng;
				exc_pos = numTS - (int)m_Op_Exc->Volt_delay[n];
				exc_pos *= (exc_pos>0);
				exc_pos %= p;
				exc_pos *= (exc_pos<(int)length);
				ny = m_Op_Exc->Volt_dir[n];
				pos[0]=m_Op_Exc->Volt_index[0][n];
				pos[1]=m_Op_Exc->Volt_index[1][n];
				pos[2]=m_Op_Exc->Volt_index[2][n];
				eng_sse->Engine_sse::SetVolt(ny,pos, eng_sse->Engine_sse::GetVolt(ny,pos) + m_Op_Exc->Volt_amp[n]*exc_volt[exc_pos]);
			}
			break;
		}
	default:
		{
			for (unsigned int n=0; n<m_Op_Exc->Volt_Count; ++n)
			{
				exc_pos = numTS - (int)m_Op_Exc->Volt_delay[n];
				exc_pos *= (exc_pos>0);
				exc_pos %= p;
				exc_pos *= (exc_pos<(int)length);
				ny = m_Op_Exc->Volt_dir[n];
				pos[0]=m_Op_Exc->Volt_index[0][n];
				pos[1]=m_Op_Exc->Volt_index[1][n];
				pos[2]=m_Op_Exc->Volt_index[2][n];
				m_Eng->SetVolt(ny,pos, m_Eng->GetVolt(ny,pos) + m_Op_Exc->Volt_amp[n]*exc_volt[exc_pos]);
			}
			break;
		}
	}
}

void Engine_Ext_Excitation::Apply2Current()
{
	//soft current excitation here (H-field excite)

	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	int numTS = m_Eng->GetNumberOfTimesteps();
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_curr =  m_Op_Exc->m_Exc->GetCurrentSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	//switch for different engine types to access faster inline engine functions
	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int n=0; n<m_Op_Exc->Curr_Count; ++n)
			{
				exc_pos = numTS - (int)m_Op_Exc->Curr_delay[n];
				exc_pos *= (exc_pos>0);
				exc_pos %= p;
				exc_pos *= (exc_pos<(int)length);
				ny = m_Op_Exc->Curr_dir[n];
				pos[0]=m_Op_Exc->Curr_index[0][n];
				pos[1]=m_Op_Exc->Curr_index[1][n];
				pos[2]=m_Op_Exc->Curr_index[2][n];
				m_Eng->Engine::SetCurr(ny,pos, m_Eng->Engine::GetCurr(ny,pos) + m_Op_Exc->Curr_amp[n]*exc_curr[exc_pos]);
			}
			break;
		}
	case Engine::SSE:
		{
			for (unsigned int n=0; n<m_Op_Exc->Curr_Count; ++n)
			{
				Engine_sse* eng_sse = (Engine_sse*) m_Eng;
				exc_pos = numTS - (int)m_Op_Exc->Curr_delay[n];
				exc_pos *= (exc_pos>0);
				exc_pos %= p;
				exc_pos *= (exc_pos<(int)length);
				ny = m_Op_Exc->Curr_dir[n];
				pos[0]=m_Op_Exc->Curr_index[0][n];
				pos[1]=m_Op_Exc->Curr_index[1][n];
				pos[2]=m_Op_Exc->Curr_index[2][n];
				eng_sse->Engine_sse::SetCurr(ny,pos, eng_sse->Engine_sse::GetCurr(ny,pos) + m_Op_Exc->Curr_amp[n]*exc_curr[exc_pos]);
			}
			break;
		}
	default:
		{
			for (unsigned int n=0; n<m_Op_Exc->Curr_Count; ++n)
			{
				exc_pos = numTS - (int)m_Op_Exc->Curr_delay[n];
				exc_pos *= (exc_pos>0);
				exc_pos %= p;
				exc_pos *= (exc_pos<(int)length);
				ny = m_Op_Exc->Curr_dir[n];
				pos[0]=m_Op_Exc->Curr_index[0][n];
				pos[1]=m_Op_Exc->Curr_index[1][n];
				pos[2]=m_Op_Exc->Curr_index[2][n];
				m_Eng->SetCurr(ny,pos, m_Eng->GetCurr(ny,pos) + m_Op_Exc->Curr_amp[n]*exc_curr[exc_pos]);
			}
			break;
		}
	}
}

// Whether the excited cell (ext_x, ext_y, ext_z) is inside the
// tile that is currently being processed.
bool
Engine_Ext_Excitation::InsideTile(
	unsigned int start[3], unsigned int stop[3],
	unsigned int ext_x, unsigned int ext_y, unsigned int ext_z
)
{
	if (ext_x < start[0] || ext_x > stop[0])
		return false;
	else if (ext_y < start[1] || ext_y > stop[1])
		return false;
	else if (ext_z < start[2] || ext_z > stop[2])
		return false;
	else
		return true;
}

void Engine_Ext_Excitation::Apply2Voltages(int timestep, unsigned int start[3], unsigned int stop[3])
{
	//soft voltage excitation here (E-field excite)
	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	int numTS = timestep;
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_volt =  m_Op_Exc->m_Exc->GetVoltageSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	if (m_Eng->GetType() != Engine::SSE)
	{
		std::cerr << "Engine_Ext_Excitation: engine unsupported for tiling!" << std::endl;
		std::exit(1);
	}

	Engine_sse* eng_sse = (Engine_sse*) m_Eng;

	for (unsigned int n=0; n<m_Op_Exc->Volt_Count; ++n)
	{
		pos[0]=m_Op_Exc->Volt_index[0][n];
		pos[1]=m_Op_Exc->Volt_index[1][n];
		pos[2]=m_Op_Exc->Volt_index[2][n];
		if (!InsideTile(start, stop, pos[0], pos[1], pos[2])) {
			continue;
		}

		exc_pos = numTS - (int)m_Op_Exc->Volt_delay[n];
		exc_pos *= (exc_pos>0);
		exc_pos %= p;
		exc_pos *= (exc_pos<(int)length);
		ny = m_Op_Exc->Volt_dir[n];
		eng_sse->Engine_sse::SetVolt(ny,pos, eng_sse->Engine_sse::GetVolt(ny,pos) + m_Op_Exc->Volt_amp[n]*exc_volt[exc_pos]);
	}
}

void Engine_Ext_Excitation::Apply2Current(int timestep, unsigned int start[3], unsigned int stop[3])
{
	//soft current excitation here (H-field excite)

	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	int numTS = timestep;
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_curr =  m_Op_Exc->m_Exc->GetCurrentSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	if (m_Eng->GetType() != Engine::SSE)
	{
		std::cerr << "Engine_Ext_Excitation: engine unsupported for tiling!" << std::endl;
		std::exit(1);
	}

	Engine_sse* eng_sse = (Engine_sse*) m_Eng;

	for (unsigned int n=0; n<m_Op_Exc->Curr_Count; ++n)
	{
		pos[0]=m_Op_Exc->Curr_index[0][n];
		pos[1]=m_Op_Exc->Curr_index[1][n];
		pos[2]=m_Op_Exc->Curr_index[2][n];
		if (!InsideTile(start, stop, pos[0], pos[1], pos[2])) {
			continue;
		}

		exc_pos = numTS - (int)m_Op_Exc->Curr_delay[n];
		exc_pos *= (exc_pos>0);
		exc_pos %= p;
		exc_pos *= (exc_pos<(int)length);
		ny = m_Op_Exc->Curr_dir[n];
		eng_sse->Engine_sse::SetCurr(ny,pos, eng_sse->Engine_sse::GetCurr(ny,pos) + m_Op_Exc->Curr_amp[n]*exc_curr[exc_pos]);
	}
}
