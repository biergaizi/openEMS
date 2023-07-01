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

#include "operator_tiling.h"
#include "engine_tiling.h"
#include "tools/useful.h"

Operator_Tiling* Operator_Tiling::New(unsigned int numThreads)
{
	cout << "Create FDTD operator (compressed SSE + multi-threading + spatial/temporal tiling)" << endl;
	cerr << "Warning! Tiling engine is highly experimental and not validated!" << endl;
	cerr << "Make sure to compare your results with the upstream openEMS for mission-critical simulations!" << endl;
	Operator_Tiling* op = new Operator_Tiling();
	op->setNumThreads(numThreads);
	op->Init();
	return op;
}

Engine* Operator_Tiling::CreateEngine()
{
	m_Engine = Engine_Tiling::New(this, m_orig_numThreads);
	return m_Engine;
}
