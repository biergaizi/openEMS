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

#include "engine_ext_dispersive.h"
#include "operator_ext_dispersive.h"
#include "FDTD/engine_sse.h"

Engine_Ext_Dispersive::Engine_Ext_Dispersive(Operator_Ext_Dispersive* op_ext_disp) : Engine_Extension(op_ext_disp)
{
	//m_TilingSupported = true;

	m_Op_Ext_Disp = op_ext_disp;
	int order = m_Op_Ext_Disp->m_Order;
	curr_ADE = new FDTD_FLOAT**[order];
	volt_ADE = new FDTD_FLOAT**[order];
	for (int o=0;o<order;++o)
	{
		curr_ADE[o] = new FDTD_FLOAT*[3];
		volt_ADE[o] = new FDTD_FLOAT*[3];
		for (int n=0; n<3; ++n)
		{
			if (m_Op_Ext_Disp->m_curr_ADE_On[o]==true)
			{
				curr_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count[o]; ++i)
					curr_ADE[o][n][i]=0.0;
			}
			else
				curr_ADE[o][n] = NULL;
			if (m_Op_Ext_Disp->m_volt_ADE_On[o]==true)
			{
				volt_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count[o]; ++i)
					volt_ADE[o][n][i]=0.0;
			}
			else
				volt_ADE[o][n] = NULL;
		}
	}
}

Engine_Ext_Dispersive::~Engine_Ext_Dispersive()
{
	if (curr_ADE==NULL && volt_ADE==NULL)
		return;

	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		for (int n=0; n<3; ++n)
		{
			delete[] curr_ADE[o][n];
			delete[] volt_ADE[o][n];
		}
		delete[] curr_ADE[o];
		delete[] volt_ADE[o];
	}
	delete[] curr_ADE;
	curr_ADE=NULL;

	delete[] volt_ADE;
	volt_ADE=NULL;

	for (auto& it : m_volt_map)
	{
		auto& vec = it.second;
		vec.clear();
	}

	for (auto& it : m_volt_map)
	{
		auto& vec = it.second;
		vec.clear();
	}
}

void Engine_Ext_Dispersive::InitializeTiling(std::vector<Range3D> tiles)
{
	int start[3];
	int end[3];

	std::cerr << "Engine_Ext_Dispersive::InitializeTiling" << std::endl;

	for (auto& tile : tiles)
	{
		for (int o = 0; o < m_Order; ++o)
		{
			if (m_Op_Ext_Disp->m_volt_ADE_On[o]==false)
				continue;

			unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];
			int* start = tile.voltageStart;
			int* end = tile.voltageStop;

			for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
			{
				if (InsideTile(start, end, pos[0][i], pos[1][i], pos[2][i]))
				{
					m_volt_map[GetTileKey(o, start, end)].push_back(i);
				}
			}
		}

		for (int o = 0; o < m_Order; ++o)
		{
			if (m_Op_Ext_Disp->m_curr_ADE_On[o]==false)
				continue;

			unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];
			int* start = tile.currentStart;
			int* end = tile.currentStop;

			for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
			{
				if (InsideTile(start, end, pos[0][i], pos[1][i], pos[2][i]))
				{
					m_curr_map[GetTileKey(o, start, end)].push_back(i);
				}
			}
		}
	}

	std::cerr << "Engine_Ext_Dispersive::InitializeTiling done" << std::endl;
}

// Whether the ADE cell (ade_x, ade_y, ade_z) is inside the
// tile that is currently being processed.
bool Engine_Ext_Dispersive::InsideTile(
	int start[3], int end[3],
	int ade_x, int ade_y, int ade_z
)
{
	int retval;

	if (ade_x < start[0] || ade_x > end[0])
		retval = false;
	else if (ade_y < start[1] || ade_y > end[1])
		retval = false;
	else if (ade_z < start[2] || ade_z > end[2])
		retval = false;
	else
		retval = true;

#if 0
	if (!retval)
	{
		fprintf(stderr, "Dispersive: cell rejected.\n");
	}
#endif
	return retval;
}

void Engine_Ext_Dispersive::Apply2Voltages(int threadID, int timestep, int start[3], int end[3])
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_volt_ADE_On[o]==false) continue;

		auto vec = m_volt_map[GetTileKey(o, start, end)];
		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		//switch for different engine types to access faster inline engine functions
		switch (m_Eng->GetType())
		{
		case Engine::BASIC:
		{
			for (auto& i : vec)
			{
				m_Eng->Engine::SetVolt(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][0][i]);
				m_Eng->Engine::SetVolt(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][1][i]);
				m_Eng->Engine::SetVolt(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][2][i]);
			}
			break;
		}
		case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*)m_Eng;
			for (auto& i : vec)
			{
				eng_sse->Engine_sse::SetVolt(0,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][0][i]);
				eng_sse->Engine_sse::SetVolt(1,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][1][i]);
				eng_sse->Engine_sse::SetVolt(2,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][2][i]);
			}
			break;
		}
		default:
			for (auto& i : vec)
			{
				m_Eng->SetVolt(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][0][i]);
				m_Eng->SetVolt(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][1][i]);
				m_Eng->SetVolt(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][2][i]);
			}
			break;
		}
	}
}

void Engine_Ext_Dispersive::Apply2Current(int threadID, int timestep, int start[3], int end[3])
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_curr_ADE_On[o]==false) continue;

		auto vec = m_curr_map[GetTileKey(o, start, end)];
		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		//switch for different engine types to access faster inline engine functions
		switch (m_Eng->GetType())
		{
		case Engine::BASIC:
		{
			for (auto& i : vec)
			{
				m_Eng->Engine::SetCurr(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][0][i]);
				m_Eng->Engine::SetCurr(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][1][i]);
				m_Eng->Engine::SetCurr(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][2][i]);
			}
			break;
		}
		case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*)m_Eng;
			for (auto& i : vec)
			{
				eng_sse->Engine_sse::SetCurr(0,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][0][i]);
				eng_sse->Engine_sse::SetCurr(1,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][1][i]);
				eng_sse->Engine_sse::SetCurr(2,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][2][i]);
			}
			break;
		}
		default:
			for (auto& i : vec)
			{
				m_Eng->SetCurr(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][0][i]);
				m_Eng->SetCurr(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][1][i]);
				m_Eng->SetCurr(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][2][i]);
			}
			break;
		}
	}
}
