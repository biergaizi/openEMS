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

#ifndef ENGINE_EXT_DISPERSIVE_H
#define ENGINE_EXT_DISPERSIVE_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "tools/TileMap.h"

class Operator_Ext_Dispersive;

class Engine_Ext_Dispersive : public Engine_Extension
{
public:
	Engine_Ext_Dispersive(Operator_Ext_Dispersive* op_ext_disp);
	virtual ~Engine_Ext_Dispersive();

	virtual void Apply2Voltages();
	virtual void Apply2Current();

	virtual void InitializeTiling(std::vector<Range3D> tiles);
	virtual void Apply2Voltages(int timestep, unsigned int start[3], unsigned int stop[3]);
	virtual void Apply2Current(int timestep, unsigned int start[3], unsigned int stop[3]);

protected:
	Operator_Ext_Dispersive* m_Op_Ext_Disp;

	//! Dispersive order
	int m_Order;

	//! ADE currents
	// Array setup: curr_ADE[N_order][direction][mesh_pos]
	FDTD_FLOAT ***curr_ADE;

	//! ADE voltages
	// Array setup: volt_ADE[N_order][direction][mesh_pos]
	FDTD_FLOAT ***volt_ADE;

	bool InsideTile(
		unsigned int start[3],
		unsigned int stop[3],
		unsigned int ade_x, unsigned int ade_y, unsigned int ade_z
	);

private:
	TileMap m_volt_map;
	TileMap m_curr_map;
};

#endif // ENGINE_EXT_DISPERSIVE_H
