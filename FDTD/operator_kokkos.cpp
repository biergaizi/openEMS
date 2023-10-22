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

#include "engine_kokkos.h"
#include "operator_kokkos.h"

Operator_Kokkos* Operator_Kokkos::New()
{
	cout << "Create FDTD operator (Kokkos)" << endl;
	Kokkos::initialize();

	Operator_Kokkos* op = new Operator_Kokkos();
	op->Init();
	return op;
}

Operator_Kokkos::Operator_Kokkos() : Operator()
{
	grid_ptr = NULL;
	vv_ptr = NULL;
	vi_ptr = NULL;
	iv_ptr = NULL;
	ii_ptr = NULL;
}

Operator_Kokkos::~Operator_Kokkos()
{
	Delete();
	Kokkos::finalize();
}

Engine* Operator_Kokkos::CreateEngine()
{
	//! create a special kokkos-engine
	m_Engine = Engine_Kokkos::New(this);
	return m_Engine;
}

void Operator_Kokkos::Init()
{
	Operator::Init();

	grid_ptr = NULL;
	vv_ptr = NULL;
	vi_ptr = NULL;
	ii_ptr = NULL;
	iv_ptr = NULL;
}

void Operator_Kokkos::Delete()
{
	delete grid_ptr;
	delete vv_ptr;
	delete vi_ptr;
	delete ii_ptr;
	delete iv_ptr;

	grid_ptr = NULL;
	vv_ptr = NULL;
	vi_ptr = NULL;
	ii_ptr = NULL;
	iv_ptr = NULL;
}

void Operator_Kokkos::Reset()
{
	Delete();
	Operator::Reset();
}


void Operator_Kokkos::InitOperator()
{
	delete grid_ptr;
	delete vv_ptr;
	delete vi_ptr;
	delete ii_ptr;
	delete iv_ptr;

	grid_ptr = new KokkosGrid(
		numLines[0], numLines[1], numLines[2],
		22, 22, 22
	);

	KokkosGrid& grid = *grid_ptr;
	vv_ptr = new KokkosGlobalArray<float>("vv", grid);
	vi_ptr = new KokkosGlobalArray<float>("vi", grid);
	ii_ptr = new KokkosGlobalArray<float>("ii", grid);
	iv_ptr = new KokkosGlobalArray<float>("iv", grid);

	fprintf(stderr, "rounding %dx%dx%d grid to %dx%dx%d (memory) and %dx%dx%d (computing)\n",
			grid.m_grid_unround_i_size, grid.m_grid_unround_j_size, grid.m_grid_unround_k_size,
			grid.m_grid_i_size, grid.m_grid_j_size, grid.m_grid_k_size,
			grid.m_grid_loadstore_i_size, grid.m_grid_loadstore_j_size, grid.m_grid_loadstore_k_size);
	fprintf(stderr, "tile size %dx%dx%d\n", grid.m_tile_i_size, grid.m_tile_j_size, grid.m_tile_k_size);
	fprintf(stderr, "sparse tile size %dx%dx%d\n",
			grid.m_sparse_tile_i_size, grid.m_sparse_tile_j_size, grid.m_sparse_tile_k_size);
	fprintf(stderr, "memory overhead = %.1f\n", (double) grid.m_grid_size / grid.m_grid_unround_size * 100.0 - 100.0);
	fprintf(stderr, "computational overhead %.1f\n", (double) grid.m_grid_loadstore_size / grid.m_grid_unround_size * 100.0 - 100.0);
	fprintf(stderr, "%d subtiles in a sparse tile\n", grid.m_sparse_subtile_num);
}
