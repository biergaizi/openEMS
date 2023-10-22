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

#ifndef SSE_CORRECT_DENORMALS
#include <xmmintrin.h>
#endif

#include "engine_kokkos.h"

// Experiments showed that the swizzled memory order no longer improves
// speed when a tile already fits in the L2 cache.
//
// TODO: remove loop swizzling or implement loop swizzling on GPU only.
#define swizzled_to_linear_lut(subtile_id)	(subtile_id)
#define linear_to_swizzled_lut(subtile_id)	(subtile_id)

using member_type = Kokkos::TeamPolicy<>::member_type;

//! \brief construct an Engine_Kokkos instance
//! it's the responsibility of the caller to free the returned pointer
Engine_Kokkos* Engine_Kokkos::New(const Operator_Kokkos* op)
{
	cout << "Create FDTD engine (Kokkos)" << endl;
	Engine_Kokkos* e = new Engine_Kokkos(op);
	e->Init();
	return e;
}

Engine_Kokkos::Engine_Kokkos(const Operator_Kokkos* op) : Engine(op)
{
	m_type = KOKKOS;
	Op = op;
	volt_ptr = NULL;
	curr_ptr = NULL;

	// speed up the calculation of denormal floating point values (flush-to-zero)
#ifndef KOKKOS_CORRECT_DENORMALS
	unsigned int oldMXCSR = _mm_getcsr(); //read the old MXCSR setting
	unsigned int newMXCSR = oldMXCSR | 0x8040; // set DAZ and FZ bits
	_mm_setcsr( newMXCSR ); //write the new MXCSR setting to the MXCSR
#endif
}

Engine_Kokkos::~Engine_Kokkos()
{
	//_mm_setcsr( oldMXCSR ); // restore old setting
	Reset();
}

void Engine_Kokkos::Init()
{
	Engine::Init();

	delete volt_ptr;
	volt_ptr = NULL; // not used

	delete curr_ptr;
	curr_ptr = NULL; // not used

	volt_ptr = new KokkosGlobalArray<float>("volt", *Op->grid_ptr);
	curr_ptr = new KokkosGlobalArray<float>("curr", *Op->grid_ptr);
}

void Engine_Kokkos::Reset()
{
	Engine::Reset();

	delete volt_ptr;
	volt_ptr = NULL; // not used

	delete curr_ptr;
	curr_ptr = NULL; // not used
}

// Profiling shows a 10%-20% speedup when forced inlining is used.
KOKKOS_INLINE_FUNCTION
__attribute__((always_inline))
static void get_neighbor_current_subtiles(
	const KokkosGrid& grid,
	const uint32_t tile_type,
	const uint32_t tile_id,
	const uint32_t linear_subtile_id,
	const KokkosGlobalArray<float>& curr_g,
	const KokkosLocalTile<float>& curr_t,
	const KokkosSubtile<float>& curr_s,
	KokkosSubtile<float>& curr_s_pi_cj_ck,
	KokkosSubtile<float>& curr_s_ci_pj_ck,
	KokkosSubtile<float>& curr_s_ci_cj_pk
)
{
	/*
	 * Using the logical linear_subtile_id, we can then find our position
	 * in the global coordinates system.
	 */
	uint32_t ti, tj, tk;
	grid.subtile_coords_to_tile(linear_subtile_id, 0, 0, 0, tile_type, ti, tj, tk);

	uint32_t gi, gj, gk;
	grid.tile_coords_to_global(tile_id, ti, tj, tk, gi, gj, gk);

	/* 
	 * Now find our three nearest neighbour subtiles.
	 * (c = center, p = previous).
	 */
	if (gi == 0) {
		/* 
		 * We are at the true boundary of the simulation box, the dependent
		 * neighbor subtiles are just the current tile itself due to truncated
		 * box.
		 */
		curr_s_pi_cj_ck = curr_s;
	}
	else if (ti == 0) {
		/* 
		 * We are at the boundary of the tile, we should get dependencies
		 * from other tiles.
		 */
		uint32_t tile_type, tile_pi_cj_ck_id, ti, tj, tk;
		grid.global_coords_to_tile(gi - 1, gj, gk, tile_type, tile_pi_cj_ck_id, ti, tj, tk);
		const KokkosTile<float>& tile_pi_cj_ck = curr_g.get_tile(tile_pi_cj_ck_id);

		uint32_t subtile_pi_cj_ck_id, subtile_pi_cj_ck_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_pi_cj_ck_linear_id, si, sj, sk);
		subtile_pi_cj_ck_id = linear_to_swizzled_lut(subtile_pi_cj_ck_linear_id);

		curr_s_pi_cj_ck = tile_pi_cj_ck.get_subtile(subtile_pi_cj_ck_id);
	}
	else {
		/* 
		 * We are at the boundary of the subtile, we should get dependencies
		 * from other subtiles.
		 */
		uint32_t subtile_pi_cj_ck_linear_id, subtile_pi_cj_ck_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti - 1, tj, tk, subtile_pi_cj_ck_linear_id, si, sj, sk);
		subtile_pi_cj_ck_id = linear_to_swizzled_lut(subtile_pi_cj_ck_linear_id);

		curr_s_pi_cj_ck = curr_t.get_subtile(subtile_pi_cj_ck_id);
	}

	if (gj == 0) {
		curr_s_ci_pj_ck = curr_s;
	}
	else if (tj == 0) {
		uint32_t tile_type, tile_ci_pj_ck_id, ti, tj, tk;
		grid.global_coords_to_tile(gi, gj - 1, gk, tile_type, tile_ci_pj_ck_id, ti, tj, tk);
		const KokkosTile<float>& tile_ci_pj_ck = curr_g.get_tile(tile_ci_pj_ck_id);

		uint32_t subtile_ci_pj_ck_id, subtile_ci_pj_ck_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_ci_pj_ck_linear_id, si, sj, sk);
		subtile_ci_pj_ck_id = linear_to_swizzled_lut(subtile_ci_pj_ck_linear_id);

		curr_s_ci_pj_ck = tile_ci_pj_ck.get_subtile(subtile_ci_pj_ck_id);
	}
	else {
		uint32_t subtile_ci_pj_ck_id, subtile_ci_pj_ck_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj - 1, tk, subtile_ci_pj_ck_linear_id, si, sj, sk);
		subtile_ci_pj_ck_id = linear_to_swizzled_lut(subtile_ci_pj_ck_linear_id);

		curr_s_ci_pj_ck = curr_t.get_subtile(subtile_ci_pj_ck_id);
	}

	if (gk == 0) {
		curr_s_ci_cj_pk = curr_s;
	}
	else if (tk == 0) {
		uint32_t tile_type, tile_ci_cj_pk_id, ti, tj, tk;
		grid.global_coords_to_tile(gi, gj, gk - 1, tile_type, tile_ci_cj_pk_id, ti, tj, tk);
		const KokkosTile<float>& tile_ci_cj_pk = curr_g.get_tile(tile_ci_cj_pk_id);

		uint32_t subtile_ci_cj_pk_id, subtile_ci_cj_pk_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_ci_cj_pk_linear_id, si, sj, sk);
		subtile_ci_cj_pk_id = linear_to_swizzled_lut(subtile_ci_cj_pk_linear_id);

		curr_s_ci_cj_pk = tile_ci_cj_pk.get_subtile(subtile_ci_cj_pk_id);
	}
	else {
		uint32_t subtile_ci_cj_pk_id, subtile_ci_cj_pk_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk - 1, subtile_ci_cj_pk_linear_id, si, sj, sk);
		subtile_ci_cj_pk_id = linear_to_swizzled_lut(subtile_ci_cj_pk_linear_id);

		curr_s_ci_cj_pk = curr_t.get_subtile(subtile_ci_cj_pk_id);
	}
}

inline static void UpdateVoltagesKernelInnerLoop(
	KokkosSubtile<float>& volt_s,
	KokkosSubtile<float>& curr_s,
	KokkosSubtile<float>& curr_s_pi_cj_ck,
	KokkosSubtile<float>& curr_s_ci_pj_ck,
	KokkosSubtile<float>& curr_s_ci_cj_pk,
	KokkosSubtile<float>& vv_s,
	KokkosSubtile<float>& vi_s,
	uint32_t gi, uint32_t gj, uint32_t gk,
	uint32_t si, uint32_t sj, uint32_t sk
)
{
	float volt0 = volt_s(0, si, sj, sk);
	float volt1 = volt_s(1, si, sj, sk);
	float volt2 = volt_s(2, si, sj, sk);

	float vv0 = vv_s(0, si, sj, sk);
	float vv1 = vv_s(1, si, sj, sk);
	float vv2 = vv_s(2, si, sj, sk);

	float vi0 = vi_s(0, si, sj, sk);
	float vi1 = vi_s(1, si, sj, sk);
	float vi2 = vi_s(2, si, sj, sk);

	float curr0_ci_cj_ck = curr_s(0, si, sj, sk);
	float curr1_ci_cj_ck = curr_s(1, si, sj, sk);
	float curr2_ci_cj_ck = curr_s(2, si, sj, sk);

	float curr1_pi_cj_ck = curr1_ci_cj_ck;
	float curr2_pi_cj_ck = curr2_ci_cj_ck;
	if (gi == 0 && si == 0) {
		/* We are at the true boundary of the simulation box, do nothing. */
	}
	else if (si == 0) {
		/* 
		 * We are at the boundary of the subtile, we should get dependencies
		 * from other subtiles.
		 */
		curr1_pi_cj_ck = curr_s_pi_cj_ck(1, SUBTILE_I_SIZE - 1, sj, sk);
		curr2_pi_cj_ck = curr_s_pi_cj_ck(2, SUBTILE_I_SIZE - 1, sj, sk);
	}
	else {
		curr1_pi_cj_ck = curr_s(1, si - 1, sj, sk);
		curr2_pi_cj_ck = curr_s(2, si - 1, sj, sk);
	}

	float curr0_ci_pj_ck = curr0_ci_cj_ck;
	float curr2_ci_pj_ck = curr2_ci_cj_ck;
	if (gj == 0 && sj == 0) {
		/* We are at the true boundary of the simulation box, do nothing. */
	}
	else if (sj == 0) {
		curr0_ci_pj_ck = curr_s_ci_pj_ck(0, si, SUBTILE_J_SIZE - 1, sk);
		curr2_ci_pj_ck = curr_s_ci_pj_ck(2, si, SUBTILE_J_SIZE - 1, sk);
	}
	else {
		curr0_ci_pj_ck = curr_s(0, si, sj - 1, sk);
		curr2_ci_pj_ck = curr_s(2, si, sj - 1, sk);
	}

	float curr0_ci_cj_pk = curr0_ci_cj_ck;
	float curr1_ci_cj_pk = curr1_ci_cj_ck;
	if (gk == 0 && sk == 0) {
		/* We are at the true boundary of the simulation box, do nothing. */
	}
	else if (sk == 0) {
		curr0_ci_cj_pk = curr_s_ci_cj_pk(0, si, sj, SUBTILE_K_SIZE - 1);
		curr1_ci_cj_pk = curr_s_ci_cj_pk(1, si, sj, SUBTILE_K_SIZE - 1);
	}
	else {
		curr0_ci_cj_pk = curr_s(0, si, sj, sk - 1);
		curr1_ci_cj_pk = curr_s(1, si, sj, sk - 1);
	}

	volt0 *= vv0;
	volt0 +=
	    vi0 * (
		curr2_ci_cj_ck -
		curr2_ci_pj_ck -
		curr1_ci_cj_ck +
		curr1_ci_cj_pk
	    );

	//for y polarization
	volt1 *= vv1;
	volt1 +=
	    vi1 * (
		curr0_ci_cj_ck -
		curr0_ci_cj_pk -
		curr2_ci_cj_ck +
		curr2_pi_cj_ck
	    );

	//for z polarization
	volt2 *= vv2;
	volt2 +=
	    vi2 * (
		curr1_ci_cj_ck -
		curr1_pi_cj_ck -
		curr0_ci_cj_ck +
		curr0_ci_pj_ck
	    );

	volt_s(0, si, sj, sk) = volt0;
	volt_s(1, si, sj, sk) = volt1;
	volt_s(2, si, sj, sk) = volt2;
}

template <bool is_sparse_tile>
static void UpdateVoltagesKernel(
	const KokkosLocalTile<float>& __restrict__ volt_t,
	const KokkosLocalTile<float>& __restrict__ curr_t,
	const KokkosLocalTile<float>& __restrict__ vv_t,
	const KokkosLocalTile<float>& __restrict__ vi_t,
	const KokkosGlobalArray<float>& __restrict__ curr_g,
	const KokkosGrid& grid,
        uint32_t tile_id,
	uint32_t tile_type
)
{
	uint32_t ti = 0, tj = 0, tk = 0;

	uint32_t gi, gj, gk;
	grid.tile_coords_to_global(tile_id, ti, tj, tk, gi, gj, gk);

	for (uint32_t subtile_id = 0; subtile_id < volt_t.subtile_num; subtile_id++) {
		KokkosSubtile<float>& volt_s = volt_t.get_subtile(subtile_id);
		KokkosSubtile<float>& curr_s = curr_t.get_subtile(subtile_id);
		KokkosSubtile<float>& vv_s = vv_t.get_subtile(subtile_id);
		KokkosSubtile<float>& vi_s = vi_t.get_subtile(subtile_id);

		/*
		 * Now the mental gymnastics begin! Although the grid and arrays
		 * are designed with a C-style linear subtile storage order, in
		 * this kernel we override it and actually read/write/iterate
		 * subtiles in a swizzled order (subtile_id) derived from some
		 * tricks in combinatorics. Thus, to know its relationships
		 * with other subtiles logically, we need the original linear
		 * subtile id.
		 */
		uint32_t linear_subtile_id;
		linear_subtile_id = swizzled_to_linear_lut(subtile_id);

		KokkosSubtile<float> curr_s_pi_cj_ck;
		KokkosSubtile<float> curr_s_ci_pj_ck;
		KokkosSubtile<float> curr_s_ci_cj_pk;
		get_neighbor_current_subtiles(
			grid,
			tile_type, tile_id, linear_subtile_id,
			curr_g, curr_t, curr_s,
			curr_s_pi_cj_ck,
			curr_s_ci_pj_ck,
			curr_s_ci_cj_pk
		);
		
		for (uint32_t si = 0; si < SUBTILE_I_SIZE; si++) {
			for (uint32_t sj = 0; sj < SUBTILE_J_SIZE; sj++) {
				for (uint32_t sk = 0; sk < SUBTILE_K_SIZE; sk++) {
					uint32_t gi, gj, gk;
					uint32_t ti, tj, tk;

					if constexpr (is_sparse_tile) {
						grid.subtile_coords_to_regular_tile(linear_subtile_id, si, sj, sk, ti, tj, tk);
					}
					else {
						grid.subtile_coords_to_sparse_tile(linear_subtile_id, si, sj, sk, tile_type, ti, tj, tk);
					}

					grid.tile_coords_to_global(tile_id, ti, tj, tk, gi, gj, gk);
					UpdateVoltagesKernelInnerLoop(
						volt_s,
						curr_s,
						curr_s_pi_cj_ck,
						curr_s_ci_pj_ck,
						curr_s_ci_cj_pk,
						vv_s,
						vi_s,
						gi, gj, gk,
						si, sj, sk
					);
				}
			}
		}
	}
}

void Engine_Kokkos::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	KokkosGrid& grid = *Op->grid_ptr;
	KokkosGlobalArray<float>& volt = *volt_ptr;
	KokkosGlobalArray<float>& curr = *curr_ptr;
	KokkosGlobalArray<float>& vv = *Op->vv_ptr;
	KokkosGlobalArray<float>& vi = *Op->vi_ptr;

	if (startX != 0 && numX != grid.m_grid_unround_i_size) {
		std::cerr << "Partial update is unimplemented in UpdateVoltages()." << std::endl;
		abort();
	}

	Kokkos::parallel_for("UpdateVoltages",
		Kokkos::TeamPolicy<>(grid.m_tile_num, 1).
			set_scratch_size(0, Kokkos::PerTeam(grid.m_tile_size * 4 * 3 * sizeof(float))),
		KOKKOS_LAMBDA (const member_type& team_member)
		{
			const uint32_t tile_id = team_member.league_rank() * team_member.team_size() +
						 team_member.team_rank();
			const uint32_t tile_type = grid.tile_id_to_tile_type(tile_id);

			KokkosLocalTile<float> scratch_volt(grid, team_member);
			KokkosLocalTile<float> scratch_curr(grid, team_member);
			KokkosLocalTile<float> scratch_vv(grid, team_member);
			KokkosLocalTile<float> scratch_vi(grid, team_member);

			scratch_volt.load_from(tile_id, volt.get_tile(tile_id));
			scratch_curr.load_from(tile_id, curr.get_tile(tile_id));
			scratch_vv.load_from(tile_id, vv.get_tile(tile_id));
			scratch_vi.load_from(tile_id, vi.get_tile(tile_id));

			switch (tile_type) {
			case TILE_REGULAR_SUBTILE:
				UpdateVoltagesKernel<false>(
					scratch_volt, scratch_curr, scratch_vv, scratch_vi, 
					curr,
					grid, tile_id, tile_type
				);
				break;
			default:
				UpdateVoltagesKernel<true>(
					scratch_volt, scratch_curr, scratch_vv, scratch_vi, 
					curr,
					grid, tile_id, tile_type
				);
				break;
			}

			scratch_volt.save_to(tile_id, volt.get_tile(tile_id));
		}
	);

	Kokkos::fence();
}

KOKKOS_INLINE_FUNCTION
__attribute__((always_inline))
static void get_neighbor_voltage_subtiles(
	const KokkosGrid& grid,
	const uint32_t tile_type,
        const uint32_t tile_id,
	const uint32_t linear_subtile_id,
	const KokkosGlobalArray<float>& volt_g,
	const KokkosLocalTile<float>& volt_t,
	const KokkosSubtile<float>& volt_s,
	KokkosSubtile<float>& volt_s_ni_cj_ck,
	KokkosSubtile<float>& volt_s_ci_nj_ck,
	KokkosSubtile<float>& volt_s_ci_cj_nk
)
{
	/*
	 * Using the logical linear_subtile_id, we can then find our position
	 * in the global coordinates system.
	 */
	uint32_t ti, tj, tk;
	grid.subtile_coords_to_tile(linear_subtile_id, 1, 1, 1, tile_type, ti, tj, tk);

	uint32_t gi, gj, gk;
	grid.tile_coords_to_global(tile_id, ti, tj, tk, gi, gj, gk);

	/*
	 * Now find our three nearest neighbour subtiles.
	 * (c = center, n = next).
	 */
	if (gi == grid.m_grid_i_size - 1) {
		/*
		 * We are at the true boundary of the simulation box, the dependent
		 * neighbor subtiles are just the current tile itself due to truncated
		 * box.
		 */
		volt_s_ni_cj_ck = volt_s;
	}
	else if (ti == volt_t.tile_i_size - 1) {
		/*
		 * We are at the boundary of the tile, we should get dependencies
		 * from other tiles.
		 */
		uint32_t tile_type, tile_ni_cj_ck_id, ti, tj, tk;
		grid.global_coords_to_tile(gi + 1, gj, gk, tile_type, tile_ni_cj_ck_id, ti, tj, tk);
		const KokkosTile<float>& tile_ni_cj_ck = volt_g.get_tile(tile_ni_cj_ck_id);

		uint32_t subtile_ni_cj_ck_id, subtile_ni_cj_ck_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_ni_cj_ck_linear_id, si, sj, sk);
		subtile_ni_cj_ck_id = linear_to_swizzled_lut(subtile_ni_cj_ck_linear_id);

		volt_s_ni_cj_ck = tile_ni_cj_ck.get_subtile(subtile_ni_cj_ck_id);
	}
	else {
		/*
		 * We are at the boundary of the subtile, we should get dependencies
		 * from other subtiles.
		 */
		uint32_t subtile_ni_cj_ck_linear_id, subtile_ni_cj_ck_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti + 1, tj, tk, subtile_ni_cj_ck_linear_id, si, sj, sk);
		subtile_ni_cj_ck_id = linear_to_swizzled_lut(subtile_ni_cj_ck_linear_id);

		volt_s_ni_cj_ck = volt_t.get_subtile(subtile_ni_cj_ck_id);
	}

	if (gj == grid.m_grid_j_size - 1) {
		volt_s_ci_nj_ck = volt_s;
	}
	else if (tj == volt_t.tile_j_size - 1) {
		uint32_t tile_type, tile_ci_nj_ck_id, ti, tj, tk;
		grid.global_coords_to_tile(gi, gj + 1, gk, tile_type, tile_ci_nj_ck_id, ti, tj, tk);
		const KokkosTile<float>& tile_ci_nj_ck = volt_g.get_tile(tile_ci_nj_ck_id);

		uint32_t subtile_ci_nj_ck_id, subtile_ci_nj_ck_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_ci_nj_ck_linear_id, si, sj, sk);
		subtile_ci_nj_ck_id = linear_to_swizzled_lut(subtile_ci_nj_ck_linear_id);

		volt_s_ci_nj_ck = tile_ci_nj_ck.get_subtile(subtile_ci_nj_ck_id);
	}
	else {
		uint32_t subtile_ci_nj_ck_id, subtile_ci_nj_ck_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj + 1, tk, subtile_ci_nj_ck_linear_id, si, sj, sk);
		subtile_ci_nj_ck_id = linear_to_swizzled_lut(subtile_ci_nj_ck_linear_id);

		volt_s_ci_nj_ck = volt_t.get_subtile(subtile_ci_nj_ck_id);
	}

	if (gk == grid.m_grid_k_size - 1) {
		volt_s_ci_cj_nk = volt_s;
	}
	else if (tk == volt_t.tile_k_size - 1) {
		uint32_t tile_type, tile_ci_cj_nk_id, ti, tj, tk;
		grid.global_coords_to_tile(gi, gj, gk + 1, tile_type, tile_ci_cj_nk_id, ti, tj, tk);
		const KokkosTile<float>& tile_ci_cj_nk = volt_g.get_tile(tile_ci_cj_nk_id);

		uint32_t subtile_ci_cj_nk_id, subtile_ci_cj_nk_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_ci_cj_nk_linear_id, si, sj, sk);
		subtile_ci_cj_nk_id = linear_to_swizzled_lut(subtile_ci_cj_nk_linear_id);

		volt_s_ci_cj_nk = tile_ci_cj_nk.get_subtile(subtile_ci_cj_nk_id);
	}
	else {
		uint32_t subtile_ci_cj_nk_id, subtile_ci_cj_nk_linear_id;
		uint32_t si, sj, sk;
		grid.tile_coords_to_subtile(tile_type, ti, tj, tk + 1, subtile_ci_cj_nk_linear_id, si, sj, sk);
		subtile_ci_cj_nk_id = linear_to_swizzled_lut(subtile_ci_cj_nk_linear_id);

		volt_s_ci_cj_nk = volt_t.get_subtile(subtile_ci_cj_nk_id);
	}
}

inline static void UpdateCurrentsKernelInnerLoop(
	KokkosSubtile<float>& curr_s,
	KokkosSubtile<float>& volt_s,
	KokkosSubtile<float>& volt_s_ni_cj_ck,
	KokkosSubtile<float>& volt_s_ci_nj_ck,
	KokkosSubtile<float>& volt_s_ci_cj_nk,
	KokkosSubtile<float>& ii_s,
	KokkosSubtile<float>& iv_s,
	uint32_t gi, uint32_t gj, uint32_t gk,
	uint32_t si, uint32_t sj, uint32_t sk,
	const KokkosGrid& grid
)
{
	float curr0 = curr_s(0, si, sj, sk);
	float curr1 = curr_s(1, si, sj, sk);
	float curr2 = curr_s(2, si, sj, sk);

	float ii0 = ii_s(0, si, sj, sk);
	float ii1 = ii_s(1, si, sj, sk);
	float ii2 = ii_s(2, si, sj, sk);

	float iv0 = iv_s(0, si, sj, sk);
	float iv1 = iv_s(1, si, sj, sk);
	float iv2 = iv_s(2, si, sj, sk);

	float volt0_ci_cj_ck = volt_s(0, si, sj, sk);
	float volt1_ci_cj_ck = volt_s(1, si, sj, sk);
	float volt2_ci_cj_ck = volt_s(2, si, sj, sk);

	float volt1_ni_cj_ck = 0;
	float volt2_ni_cj_ck = 0;
	if (gi == grid.m_grid_unround_i_size - 1) {
		/*
		 * dependencies are outside simulation box (boundary condition),
		 * replace operators with no-op.
		 */
		ii0 = 1;
		ii1 = 1;
		ii2 = 1;
		iv0 = 0;
		iv1 = 0;
		iv2 = 0;
	}
	else if (si == 0) {
		volt1_ni_cj_ck = volt_s(1, si + 1, sj, sk);
		volt2_ni_cj_ck = volt_s(2, si + 1, sj, sk);
	}
	else {
		volt1_ni_cj_ck = volt_s_ni_cj_ck(1, 0, sj, sk);
		volt2_ni_cj_ck = volt_s_ni_cj_ck(2, 0, sj, sk);
	}

	float volt0_ci_nj_ck = 0;
	float volt2_ci_nj_ck = 0;
	if (gj == grid.m_grid_unround_j_size - 1) {
		/*
		 * dependencies are outside simulation box (boundary condition),
		 * replace operators with no-op.
		 */
		ii0 = 1;
		ii1 = 1;
		ii2 = 1;
		iv0 = 0;
		iv1 = 0;
		iv2 = 0;
	}
	else if (sj == 0) {
		volt0_ci_nj_ck = volt_s(0, si, sj + 1, sk);
		volt2_ci_nj_ck = volt_s(2, si, sj + 1, sk);
	}
	else {
		volt0_ci_nj_ck = volt_s_ci_nj_ck(0, si, 0, sk);
		volt2_ci_nj_ck = volt_s_ci_nj_ck(2, si, 0, sk);
	}

	float volt0_ci_cj_nk = 0;
	float volt1_ci_cj_nk = 0;
	if (gk == grid.m_grid_unround_k_size - 1) {
		/*
		 * dependencies are outside simulation box (boundary condition),
		 * replace operators with no-op.
		 */
		ii0 = 1;
		ii1 = 1;
		ii2 = 1;
		iv0 = 0;
		iv1 = 0;
		iv2 = 0;
	}
	else if (sk == 0) {
		volt0_ci_cj_nk = volt_s(0, si, sj, sk + 1);
		volt1_ci_cj_nk = volt_s(1, si, sj, sk + 1);
	}
	else {
		volt0_ci_cj_nk = volt_s_ci_cj_nk(0, si, sj, 0);
		volt1_ci_cj_nk = volt_s_ci_cj_nk(1, si, sj, 0);
	}

	//for x polarization
	curr0 *= ii0;
	curr0 +=
	    iv0 * (
		volt2_ci_cj_ck -
		volt2_ci_nj_ck -
		volt1_ci_cj_ck +
		volt1_ci_cj_nk
	    );

	//for y polarization
	curr1 *= ii1;
	curr1 +=
	    iv1 * (
		volt0_ci_cj_ck -
		volt0_ci_cj_nk -
		volt2_ci_cj_ck +
		volt2_ni_cj_ck
	    );

	//for z polarization
	curr2 *= ii2;
	curr2 +=
	    iv2 * (
		volt1_ci_cj_ck -
		volt1_ni_cj_ck -
		volt0_ci_cj_ck +
		volt0_ci_nj_ck
	    );

	curr_s(0, si, sj, sk) = curr0;
	curr_s(1, si, sj, sk) = curr1;
	curr_s(2, si, sj, sk) = curr2;
}

template <bool is_sparse_tile>
static void UpdateCurrentsKernel(
	KokkosLocalTile<float>& __restrict__ curr_t,
	KokkosLocalTile<float>& __restrict__ volt_t,
	KokkosLocalTile<float>& __restrict__ ii_t,
	KokkosLocalTile<float>& __restrict__ iv_t,
	const KokkosGlobalArray<float>& __restrict__ volt_g,
	const KokkosGrid& grid,
        uint32_t tile_id,
	uint32_t tile_type
)
{
	for (uint32_t subtile_id = 0; subtile_id < curr_t.subtile_num; subtile_id++) {
		KokkosSubtile<float>& curr_s = curr_t.get_subtile(subtile_id);
		KokkosSubtile<float>& volt_s = volt_t.get_subtile(subtile_id);
		KokkosSubtile<float>& ii_s = ii_t.get_subtile(subtile_id);
		KokkosSubtile<float>& iv_s = iv_t.get_subtile(subtile_id);

		/*
		 * Now the mental gymnastics begin! Although the grid and arrays
		 * are designed with a C-style linear subtile storage order, in
		 * this kernel we override it and actually read/write/iterate
		 * subtiles in a swizzled order (subtile_id) derived from some
		 * tricks in combinatorics. Thus, to know its relationships
		 * with other subtiles logically, we need the original linear
		 * subtile id.
		 */
		uint32_t linear_subtile_id;
		linear_subtile_id = swizzled_to_linear_lut(subtile_id);

		KokkosSubtile<float> volt_s_ni_cj_ck;
		KokkosSubtile<float> volt_s_ci_nj_ck;
		KokkosSubtile<float> volt_s_ci_cj_nk;
		get_neighbor_voltage_subtiles(
			grid,
			tile_type, tile_id, linear_subtile_id,
			volt_g, volt_t, volt_s,
			volt_s_ni_cj_ck,
			volt_s_ci_nj_ck,
			volt_s_ci_cj_nk
		);

		for (uint32_t si = 0; si < SUBTILE_I_SIZE; si++) {
			for (uint32_t sj = 0; sj < SUBTILE_J_SIZE; sj++) {
				for (uint32_t sk = 0; sk < SUBTILE_K_SIZE; sk++) {
					uint32_t gi, gj, gk;
					uint32_t ti, tj, tk;

					if constexpr (is_sparse_tile) {
						grid.subtile_coords_to_regular_tile(linear_subtile_id, si, sj, sk, ti, tj, tk);
					}
					else {
						grid.subtile_coords_to_sparse_tile(linear_subtile_id, si, sj, sk, tile_type, ti, tj, tk);
					}
					grid.tile_coords_to_global(tile_id, ti, tj, tk, gi, gj, gk);

					UpdateCurrentsKernelInnerLoop(
						curr_s,
						volt_s,
						volt_s_ni_cj_ck,
						volt_s_ci_nj_ck,
						volt_s_ci_cj_nk,
						ii_s,
						iv_s,
						gi, gj, gk,
						si, sj, sk,
						grid
					);
				}
			}
		}
	}
}

void Engine_Kokkos::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	KokkosGrid& grid = *Op->grid_ptr;
	KokkosGlobalArray<float>& curr = *curr_ptr;
	KokkosGlobalArray<float>& volt = *volt_ptr;
	KokkosGlobalArray<float>& ii = *Op->ii_ptr;
	KokkosGlobalArray<float>& iv = *Op->iv_ptr;

	if (startX != 0 && numX != grid.m_grid_unround_i_size) {
		std::cerr << "Partial update is unimplemented in UpdateCurrents()." << std::endl;
		abort();
	}

	Kokkos::parallel_for("UpdateCurrents",
		Kokkos::TeamPolicy<>(grid.m_tile_num, 1).
			set_scratch_size(0, Kokkos::PerTeam(grid.m_tile_size * 4 * 3 * sizeof(float))),
		KOKKOS_LAMBDA (const member_type& team_member)
		{
			const uint32_t tile_id = team_member.league_rank() * team_member.team_size() +
						 team_member.team_rank();
			const uint32_t tile_type = grid.tile_id_to_tile_type(tile_id);

			KokkosLocalTile<float> scratch_curr(grid, team_member);
			KokkosLocalTile<float> scratch_volt(grid, team_member);
			KokkosLocalTile<float> scratch_ii(grid, team_member);
			KokkosLocalTile<float> scratch_iv(grid, team_member);

			scratch_curr.load_from(tile_id, curr.get_tile(tile_id));
			scratch_volt.load_from(tile_id, volt.get_tile(tile_id));
			scratch_ii.load_from(tile_id, ii.get_tile(tile_id));
			scratch_iv.load_from(tile_id, iv.get_tile(tile_id));

			switch (tile_type) {
			case TILE_REGULAR_SUBTILE:
				UpdateCurrentsKernel<false>(
					scratch_curr, scratch_volt, scratch_ii, scratch_iv, 
					volt,
					grid, tile_id, tile_type
				);
				break;
			default:
				UpdateCurrentsKernel<true>(
					scratch_curr, scratch_volt, scratch_ii, scratch_iv, 
					volt,
					grid, tile_id, tile_type
				);
				break;
			}

			scratch_curr.save_to(tile_id, curr.get_tile(tile_id));
		}
	);

	Kokkos::fence();
}
