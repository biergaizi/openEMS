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

#include "engine_ext_upml.h"
#include "operator_ext_upml.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/useful.h"

Engine_Ext_UPML::Engine_Ext_UPML(Operator_Ext_UPML* op_ext) : Engine_Extension(op_ext)
{
	m_Op_UPML = op_ext;

	//this ABC extension should be executed first!
	m_Priority = ENG_EXT_PRIO_UPML;

	// this extension support the tiling engine
	m_TilingSupported = true;

	volt_flux_ptr = Create_Flat_N_3DArray<FDTD_FLOAT>(m_Op_UPML->m_numLines);
	curr_flux_ptr = Create_Flat_N_3DArray<FDTD_FLOAT>(m_Op_UPML->m_numLines);

	// only matters for non-tiling engines, in the tiling engine,
	// the domain partitioning is determined by the engine, not
	// the extension
	SetNumberOfThreads(1);
}

Engine_Ext_UPML::~Engine_Ext_UPML()
{
	Delete_Flat_N_3DArray(volt_flux_ptr,m_Op_UPML->m_numLines);
	volt_flux_ptr=NULL;
	Delete_Flat_N_3DArray(curr_flux_ptr,m_Op_UPML->m_numLines);
	curr_flux_ptr=NULL;
}

void Engine_Ext_UPML::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	m_numX = AssignJobs2Threads(m_Op_UPML->m_numLines[0],m_NrThreads,false);
	m_start.resize(m_NrThreads,0);
	m_start.at(0)=0;
	for (size_t n=1; n<m_numX.size(); ++n)
		m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}

// When the tiling engine is used (methods with arguments "timestep",
// "start", and "stop"), the global 3D space is divided into small
// tiles for processing. We need to know whether the current tile
// overlaps with UPML. If there's no overlap, return immediately.
// If there's overlap, return the local UPML coordinates of the
// overlapped region.
bool
Engine_Ext_UPML::ToLocalCoords(
	unsigned int start[3], unsigned int end[3],
	unsigned int localStart[3], unsigned int localEnd[3]
)
{
	unsigned int pmlEnd[3] = {
		m_Op_UPML->m_StartPos[0] + m_Op_UPML->m_numLines[0] - 1,
		m_Op_UPML->m_StartPos[1] + m_Op_UPML->m_numLines[1] - 1,
		m_Op_UPML->m_StartPos[2] + m_Op_UPML->m_numLines[2] - 1
	};

	if (!start && !end)
	{
		const static int threadID = 0;
		localStart[0] = m_start.at(0);
		localStart[1] = 0;
		localStart[2] = 0;

		localEnd[0] = m_numX.at(0) - 1;
		localEnd[1] = m_Op_UPML->m_numLines[1] - 1;
		localEnd[2] = m_Op_UPML->m_numLines[2] - 1;

		return false;
	}

	for (int i = 0; i < 3; i++)
	{
		if (start[i] <= pmlEnd[i] && end[i] >= m_Op_UPML->m_StartPos[i])
		{
			// Return the overlapped region between the current tile and UPML>
			localStart[i] = std::max((int) start[i], (int) m_Op_UPML->m_StartPos[i]);
			localEnd[i] = std::min((int) end[i], (int) pmlEnd[i]);
			localEnd[i] -= m_Op_UPML->m_StartPos[i];
			localStart[i] -= m_Op_UPML->m_StartPos[i];
		}
		else
		{
			// Tile and UPML do not overlap.
			return false;
		}
	}

	return true;
}

void Engine_Ext_UPML::DoPreVoltageUpdatesKernel(unsigned int loc_start[3], unsigned int loc_stop[3])
{
	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	Flat_N_3DArray<FDTD_FLOAT>& volt_flux = *volt_flux_ptr;
	Flat_N_3DArray<FDTD_FLOAT>& op_upml_vv = *(m_Op_UPML->vv_ptr);
	Flat_N_3DArray<FDTD_FLOAT>& op_upml_vvfo = *(m_Op_UPML->vvfo_ptr);

	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = op_upml_vv(0, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->Engine::GetVolt(0,pos)
						         - op_upml_vvfo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->Engine::SetVolt(0,pos, volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_vv(1, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->Engine::GetVolt(1,pos)
						         - op_upml_vvfo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->Engine::SetVolt(1,pos, volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_vv(2, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->Engine::GetVolt(2,pos)
						         - op_upml_vvfo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->Engine::SetVolt(2,pos, volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = op_upml_vv(0, loc_pos[0], loc_pos[1], loc_pos[2])   * eng_sse->Engine_sse::GetVolt(0,pos)
						         - op_upml_vvfo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						eng_sse->Engine_sse::SetVolt(0,pos, volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_vv(1, loc_pos[0], loc_pos[1], loc_pos[2])   * eng_sse->Engine_sse::GetVolt(1,pos)
						         - op_upml_vvfo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						eng_sse->Engine_sse::SetVolt(1,pos, volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_vv(2, loc_pos[0], loc_pos[1], loc_pos[2])   * eng_sse->Engine_sse::GetVolt(2,pos)
						         - op_upml_vvfo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						eng_sse->Engine_sse::SetVolt(2,pos, volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;
					}
				}
			}
			break;
		}
	default:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = op_upml_vv(0, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->GetVolt(0,pos)
						         - op_upml_vvfo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->SetVolt(0,pos, volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_vv(1, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->GetVolt(1,pos)
						         - op_upml_vvfo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->SetVolt(1,pos, volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_vv(2, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->GetVolt(2,pos)
						         - op_upml_vvfo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->SetVolt(2,pos, volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
						volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;
					}
				}
			}
			break;
		}
	}

}

void Engine_Ext_UPML::DoPreVoltageUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	// Calculate the start and stop lines using UPML's local
	// coordinate system. For example, if UPML exists between
	// (10, 10, 10) and (18, 18, 18) in the main grid, for UPML
	// it's (0, 0, 0) to (8, 8, 8).

	unsigned int loc_start[3] = {
		0, 0, 0
	};
	unsigned int loc_stop[3] = {
		m_numX.at(threadID) - 1,
		m_Op_UPML->m_numLines[1] - 1,
		m_Op_UPML->m_numLines[2] - 1
	};

	DoPreVoltageUpdatesKernel(loc_start, loc_stop);
}

void Engine_Ext_UPML::DoPreVoltageUpdates(int timestep, unsigned int start[3], unsigned int stop[3])
{
	unsigned int loc_start[3];
	unsigned int loc_stop[3];

	bool insidePML = ToLocalCoords(start, stop, loc_start, loc_stop);
	if (!insidePML)
	{
		// no overlap between the current tile and UPML, do nothing.
		return;
	}

	DoPreVoltageUpdatesKernel(loc_start, loc_stop);
}

void Engine_Ext_UPML::DoPostVoltageUpdatesKernel(unsigned int loc_start[3], unsigned int loc_stop[3])
{
	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	Flat_N_3DArray<FDTD_FLOAT>& volt_flux = *volt_flux_ptr;
	Flat_N_3DArray<FDTD_FLOAT>& op_upml_vvfn = *(m_Op_UPML->vvfn_ptr);

	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->Engine::GetVolt(0,pos);
						m_Eng->Engine::SetVolt(0,pos, f_help + op_upml_vvfn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->Engine::GetVolt(1,pos);
						m_Eng->Engine::SetVolt(1,pos, f_help + op_upml_vvfn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->Engine::GetVolt(2,pos);
						m_Eng->Engine::SetVolt(2,pos, f_help + op_upml_vvfn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = eng_sse->Engine_sse::GetVolt(0,pos);
						eng_sse->Engine_sse::SetVolt(0,pos, f_help + op_upml_vvfn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = eng_sse->Engine_sse::GetVolt(1,pos);
						eng_sse->Engine_sse::SetVolt(1,pos, f_help + op_upml_vvfn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = eng_sse->Engine_sse::GetVolt(2,pos);
						eng_sse->Engine_sse::SetVolt(2,pos, f_help + op_upml_vvfn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
					}
				}
			}
			break;
		}
	default:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->GetVolt(0,pos);
						m_Eng->SetVolt(0,pos, f_help + op_upml_vvfn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->GetVolt(1,pos);
						m_Eng->SetVolt(1,pos, f_help + op_upml_vvfn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->GetVolt(2,pos);
						m_Eng->SetVolt(2,pos, f_help + op_upml_vvfn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
					}
				}
			}
			break;
		}
	}

}

void Engine_Ext_UPML::DoPostVoltageUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	// Calculate the start and stop lines using UPML's local
	// coordinate system. For example, if UPML exists between
	// (10, 10, 10) and (18, 18, 18) in the main grid, for UPML
	// it's (0, 0, 0) to (8, 8, 8).

	unsigned int loc_start[3] = {
		0, 0, 0
	};
	unsigned int loc_stop[3] = {
		m_numX.at(threadID) - 1,
		m_Op_UPML->m_numLines[1] - 1,
		m_Op_UPML->m_numLines[2] - 1
	};

	DoPostVoltageUpdatesKernel(loc_start, loc_stop);
}

void Engine_Ext_UPML::DoPostVoltageUpdates(int timestep, unsigned int start[3], unsigned int stop[3])
{
	unsigned int loc_start[3];
	unsigned int loc_stop[3];

	bool insidePML = ToLocalCoords(start, stop, loc_start, loc_stop);
	if (!insidePML)
	{
		// no overlap between the current tile and UPML, do nothing.
		return;
	}

	DoPostVoltageUpdatesKernel(loc_start, loc_stop);
}

void Engine_Ext_UPML::DoPreCurrentUpdatesKernel(unsigned int loc_start[3], unsigned int loc_stop[3])
{
	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	Flat_N_3DArray<FDTD_FLOAT>& curr_flux = *curr_flux_ptr;
	Flat_N_3DArray<FDTD_FLOAT>& op_upml_ii = *(m_Op_UPML->ii_ptr);
	Flat_N_3DArray<FDTD_FLOAT>& op_upml_iifo = *(m_Op_UPML->iifo_ptr);

	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = op_upml_ii(0, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->Engine::GetCurr(0,pos)
						         - op_upml_iifo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->Engine::SetCurr(0,pos, curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_ii(1, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->Engine::GetCurr(1,pos)
						         - op_upml_iifo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->Engine::SetCurr(1,pos, curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_ii(2, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->Engine::GetCurr(2,pos)
						         - op_upml_iifo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->Engine::SetCurr(2,pos, curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = op_upml_ii(0, loc_pos[0], loc_pos[1], loc_pos[2])   * eng_sse->Engine_sse::GetCurr(0,pos)
						         - op_upml_iifo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						eng_sse->Engine_sse::SetCurr(0,pos, curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_ii(1, loc_pos[0], loc_pos[1], loc_pos[2])   * eng_sse->Engine_sse::GetCurr(1,pos)
						         - op_upml_iifo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						eng_sse->Engine_sse::SetCurr(1,pos, curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_ii(2, loc_pos[0], loc_pos[1], loc_pos[2])   * eng_sse->Engine_sse::GetCurr(2,pos)
						         - op_upml_iifo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						eng_sse->Engine_sse::SetCurr(2,pos, curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

					}
				}
			}
			break;
		}
	default:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = op_upml_ii(0, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->GetCurr(0,pos)
						         - op_upml_iifo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->SetCurr(0,pos, curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_ii(1, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->GetCurr(1,pos)
						         - op_upml_iifo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->SetCurr(1,pos, curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

						f_help = op_upml_ii(2, loc_pos[0], loc_pos[1], loc_pos[2])   * m_Eng->GetCurr(2,pos)
						         - op_upml_iifo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						m_Eng->SetCurr(2,pos, curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
						curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;
					}
				}
			}
			break;
		}
	}
}

void Engine_Ext_UPML::DoPreCurrentUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	// Calculate the start and stop lines using UPML's local
	// coordinate system. For example, if UPML exists between
	// (10, 10, 10) and (18, 18, 18) in the main grid, for UPML
	// it's (0, 0, 0) to (8, 8, 8).

	unsigned int loc_start[3] = {
		0, 0, 0
	};
	unsigned int loc_stop[3] = {
		m_numX.at(threadID) - 1,
		m_Op_UPML->m_numLines[1] - 1,
		m_Op_UPML->m_numLines[2] - 1
	};

	DoPreCurrentUpdatesKernel(loc_start, loc_stop);
}

void Engine_Ext_UPML::DoPreCurrentUpdates(int timestep, unsigned int start[3], unsigned int stop[3])
{
	unsigned int loc_start[3];
	unsigned int loc_stop[3];

	bool insidePML = ToLocalCoords(start, stop, loc_start, loc_stop);
	if (!insidePML)
	{
		// no overlap between the current tile and UPML, do nothing.
		return;
	}

	DoPreCurrentUpdatesKernel(loc_start, loc_stop);
}

void Engine_Ext_UPML::DoPostCurrentUpdatesKernel(unsigned int loc_start[3], unsigned int loc_stop[3])
{
	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	Flat_N_3DArray<FDTD_FLOAT>& curr_flux = *curr_flux_ptr;
	Flat_N_3DArray<FDTD_FLOAT>& op_upml_ii = *(m_Op_UPML->ii_ptr);
	Flat_N_3DArray<FDTD_FLOAT>& op_upml_iifn = *(m_Op_UPML->iifn_ptr);

	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->Engine::GetCurr(0,pos);
						m_Eng->Engine::SetCurr(0,pos, f_help + op_upml_iifn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->Engine::GetCurr(1,pos);
						m_Eng->Engine::SetCurr(1,pos, f_help + op_upml_iifn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->Engine::GetCurr(2,pos);
						m_Eng->Engine::SetCurr(2,pos, f_help + op_upml_iifn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = eng_sse->Engine_sse::GetCurr(0,pos);
						eng_sse->Engine_sse::SetCurr(0,pos, f_help + op_upml_iifn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = eng_sse->Engine_sse::GetCurr(1,pos);
						eng_sse->Engine_sse::SetCurr(1,pos, f_help + op_upml_iifn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = eng_sse->Engine_sse::GetCurr(2,pos);
						eng_sse->Engine_sse::SetCurr(2,pos, f_help + op_upml_iifn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
					}
				}
			}
			break;
		}
	default:
		{
			for (loc_pos[0] = loc_start[0]; loc_pos[0] <= loc_stop[0]; ++loc_pos[0])
			{
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1] = loc_start[1]; loc_pos[1] <= loc_stop[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2] = loc_start[2]; loc_pos[2] <= loc_stop[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->GetCurr(0,pos);
						m_Eng->SetCurr(0,pos, f_help + op_upml_iifn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->GetCurr(1,pos);
						m_Eng->SetCurr(1,pos, f_help + op_upml_iifn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

						f_help = curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
						curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = m_Eng->GetCurr(2,pos);
						m_Eng->SetCurr(2,pos, f_help + op_upml_iifn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
					}
				}
			}
			break;
		}
	}
}

void Engine_Ext_UPML::DoPostCurrentUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	// Calculate the start and stop lines using UPML's local
	// coordinate system. For example, if UPML exists between
	// (10, 10, 10) and (18, 18, 18) in the main grid, for UPML
	// it's (0, 0, 0) to (8, 8, 8).

	unsigned int loc_start[3] = {
		0, 0, 0
	};
	unsigned int loc_stop[3] = {
		m_numX.at(threadID) - 1,
		m_Op_UPML->m_numLines[1] - 1,
		m_Op_UPML->m_numLines[2] - 1
	};

	DoPostCurrentUpdatesKernel(loc_start, loc_stop);
}

void Engine_Ext_UPML::DoPostCurrentUpdates(int timestep, unsigned int start[3], unsigned int stop[3])
{
	unsigned int loc_start[3];
	unsigned int loc_stop[3];

	bool insidePML = ToLocalCoords(start, stop, loc_start, loc_stop);
	if (!insidePML)
	{
		// no overlap between the current tile and UPML, do nothing.
		return;
	}

	DoPostCurrentUpdatesKernel(loc_start, loc_stop);
}
