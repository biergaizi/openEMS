/*
 * Copyright (C) 2023 Yifeng Li (tomli@tomli.me)
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
#include <Kokkos_Core.hpp>

/*
 * A subtile is the basic unit of tiling in this engine's tiling scheme.
 * It has a fixed size of 2x2x2. A tile is a combinations of multiple
 * subtiles. This arrangement serves two purposes - (1) During operator
 * compression, duplicated operators are eliminated at the subtile
 * granularity. (2) If a tile is too large to fit into on-chip memory
 * but the use of that size is desirable for other reasons, it can be
 * processed in multiple batches, with each batch fits in cache. Both
 * uses are currently unimplemented.
 */
#define SUBTILE_I_SIZE		2
#define SUBTILE_J_SIZE		2
#define SUBTILE_K_SIZE		2
#define SUBTILE_SIZE		(SUBTILE_I_SIZE * SUBTILE_J_SIZE * SUBTILE_K_SIZE)

template <typename T, uint32_t n_max=3>
struct KokkosSubtile
{
	T array[n_max * SUBTILE_SIZE];

	/*
	 * Convert a cell from the subtile coordinates from the perspective
	 * within a 2x2x2 subtile (si, sj, sk) to its linear address from
	 * the perspective of a 1D array.
	 *
	 * Don't forget to update grid.subtile_coords_to_linear() if this
	 * has changed.
	 */
	inline T& operator() (
		uint32_t n,
		uint32_t si,
		uint32_t sj,
		uint32_t sk
	)
	{
		const uint32_t n_size = 3;

		const uint32_t i_stride = SUBTILE_J_SIZE * SUBTILE_K_SIZE * n_size;
		const uint32_t j_stride = SUBTILE_K_SIZE * n_size;
		const uint32_t k_stride = n_size;

		return array[si * i_stride + sj * j_stride + sk * k_stride + n];
	}
};

/*
 * A tile is the basic unit of computation. During simulation, tiles
 * are loaded one by one from memory, and all processing is done
 * on a tile. By default, a tile has a size of 32x32x32 (and contains
 * 4096 subtiles), but experiments shows best performance when a tile
 * fits in a CPU's L2 cache, which usually has a smaller size. Thus,
 * the size can be adjusted on runtime.
 */
#define TILE_I_DEFAULT_SIZE	32
#define TILE_J_DEFAULT_SIZE	32
#define TILE_K_DEFAULT_SIZE	32

enum {
	TILE_REGULAR_SUBTILE = 0,
	TILE_SPARSE_SUBTILE = 1,

	/*
	 * The I, J, K dimensions can be a short dimension or a full dimension
	 * on its own, so there are 7 possible combinations in total. What a
	 * curse of dimensionality...
	 */
	TILE_SPARSE_SUBTILE_I = 2,
	TILE_SPARSE_SUBTILE_J = 4,
	TILE_SPARSE_SUBTILE_K = 8,
};

/*
 * Grid is the manager of size and coordinates transformation at all
 * levels, including the entire simulation box, a tile, or a subtile.
 * Before a simulation, a grid is first created with a specific domain
 * size and tile size. Then, all the arrays are created using the
 * information from the grid. If a simulation domain is not a multiple
 * of the tile size, it should also do the necessary roundings (not
 * implemented).
 */
struct KokkosGrid
{
	KokkosGrid(
		uint32_t grid_i_size,
		uint32_t grid_j_size,
		uint32_t grid_k_size,
		uint32_t tile_i_size=TILE_I_DEFAULT_SIZE,
		uint32_t tile_j_size=TILE_J_DEFAULT_SIZE,
		uint32_t tile_k_size=TILE_K_DEFAULT_SIZE
	)
	{
		rounded = false;
		m_grid_unround_i_size = grid_i_size;
		m_grid_unround_j_size = grid_j_size;
		m_grid_unround_k_size = grid_k_size;
		m_grid_unround_size = grid_i_size * grid_j_size * grid_k_size;

		/*
		 * Each dimension of a tile must be padding to an even
		 * number, because the smallest subtile is 2x2x2. This
		 * padding overhead always exists.
		 */
		bool tile_rounded = false;
		if (tile_i_size % 2 != 0) {
			tile_i_size++;
			tile_rounded = true;
		}
		if (tile_j_size % 2 != 0) {
			tile_j_size++;
			tile_rounded = true;
		}
		if (tile_k_size % 2 != 0) {
			tile_k_size++;
			tile_rounded = true;
		}
		if (tile_rounded) {
			fprintf(stderr, "warning: rounded tiling size should be %dx%dx%d\n",
					tile_i_size, tile_i_size, tile_k_size);
			exit(1);
		}

		/*
		 * If the simulation box is not an integer multiple of a tile's
		 * size, we pad the simulation box's size to an integer multiple
		 * of a tile's size. But note that they're created using create
		 * special "sparse tile" at the edges of the simulation box. Memory
		 * allocation and global-to-tile coordinate calculation still act
		 * as if these sparse tiles are full tiles. But actual reads/writes
		 * to these sparse tiles are much smaller than a full tile by using
		 * a different formula for tile-to-subtile coordinate calculation.
		 * Instead of paying the full pad-to-tile-size overhead up to 30%,
		 * they only have a pad-to-subtile-size overhead around 1%.
		 */
		if (grid_i_size % tile_i_size != 0) {
			grid_i_size += tile_i_size - grid_i_size % tile_i_size;
			rounded = true;
		}
		if (grid_j_size % tile_j_size != 0) {
			grid_j_size += tile_j_size - grid_j_size % tile_j_size;
			rounded = true;
		}
		if (grid_k_size % tile_k_size != 0) {
			grid_k_size += tile_k_size - grid_k_size % tile_k_size;
			rounded = true;
		}

		/*
		 * TODO: remove all m_ prefixes from member variables that
		 * are meant to be publically accessible.
		 */

		/* how large is the simulation box, i.e. how many cells? */
		m_grid_i_size = grid_i_size;
		m_grid_j_size = grid_j_size;
		m_grid_k_size = grid_k_size;
		m_grid_size = grid_i_size * grid_j_size * grid_k_size;

		/* how large is a standard tile, i.e. how many cells?? */
		m_tile_i_size = tile_i_size;
		m_tile_j_size = tile_j_size;
		m_tile_k_size = tile_k_size;
		m_tile_size = tile_i_size * tile_j_size * tile_k_size;

		/* how many standard tiles are there in the simulation box? */
		m_tile_i_num = grid_i_size / tile_i_size;
		m_tile_j_num = grid_j_size / tile_j_size;
		m_tile_k_num = grid_k_size / tile_k_size;
		m_tile_num = m_grid_size / m_tile_size;

		/* how many subtiles are there in a standard tile? */
		m_subtile_i_num = tile_i_size / SUBTILE_I_SIZE;
		m_subtile_j_num = tile_j_size / SUBTILE_J_SIZE;
		m_subtile_k_num = tile_k_size / SUBTILE_K_SIZE;
		m_subtile_num = m_tile_size / SUBTILE_SIZE;

		/* how large is a sparse tile? */
		m_sparse_tile_i_size = m_tile_i_size - (m_grid_i_size - m_grid_unround_i_size);
		if (m_sparse_tile_i_size % SUBTILE_I_SIZE != 0) {
			/*
			 * a sparse tile must be divisible by 2, since the smallest unit is still
			 * a 2x2x2 subtile.
			 */
			m_sparse_tile_i_size++;
		}

		m_sparse_tile_j_size = m_tile_j_size - (m_grid_j_size - m_grid_unround_j_size);
		if (m_sparse_tile_j_size % SUBTILE_J_SIZE != 0) {
			m_sparse_tile_j_size++;
		}

		m_sparse_tile_k_size = m_tile_k_size - (m_grid_k_size - m_grid_unround_k_size);
		if (m_sparse_tile_k_size % SUBTILE_K_SIZE != 0) {
			m_sparse_tile_k_size++;
		}

		m_sparse_tile_size = m_sparse_tile_i_size * m_sparse_tile_j_size * m_sparse_tile_k_size;

		/* how many subtiles are there in a sparse tile? */
		m_sparse_subtile_i_num = m_sparse_tile_i_size / SUBTILE_I_SIZE;
		m_sparse_subtile_j_num = m_sparse_tile_j_size / SUBTILE_J_SIZE;
		m_sparse_subtile_k_num = m_sparse_tile_k_size / SUBTILE_K_SIZE;
		m_sparse_subtile_num = m_sparse_tile_size / SUBTILE_SIZE;

		m_grid_loadstore_i_size = (m_tile_i_num - 1) * m_tile_i_size + m_sparse_tile_i_size;
		m_grid_loadstore_j_size = (m_tile_j_num - 1) * m_tile_j_size + m_sparse_tile_j_size;
		m_grid_loadstore_k_size = (m_tile_k_num - 1) * m_tile_k_size + m_sparse_tile_k_size;
		m_grid_loadstore_size = m_grid_loadstore_i_size * m_grid_unround_j_size * m_grid_loadstore_k_size;
	};

	bool rounded;
	uint32_t m_grid_unround_i_size, m_grid_unround_j_size, m_grid_unround_k_size, m_grid_unround_size;
	uint32_t m_grid_loadstore_i_size, m_grid_loadstore_j_size, m_grid_loadstore_k_size, m_grid_loadstore_size;

	uint32_t m_grid_i_size, m_grid_j_size, m_grid_k_size, m_grid_size;
	uint32_t m_tile_i_size, m_tile_j_size, m_tile_k_size, m_tile_size;
	uint32_t m_tile_i_num, m_tile_j_num, m_tile_k_num, m_tile_num;
	uint32_t m_subtile_i_num, m_subtile_j_num, m_subtile_k_num, m_subtile_num;
	uint32_t m_sparse_tile_i_size, m_sparse_tile_j_size, m_sparse_tile_k_size, m_sparse_tile_size;
	uint32_t m_sparse_subtile_i_num, m_sparse_subtile_j_num, m_sparse_subtile_k_num, m_sparse_subtile_num;

	/*
	 * Convert a cell from the global coordinates from the perspective
	 * within a global simulatino box (gi, gj, gk) to its tile coordinates
	 * from the perspective of a 32x32x32 tile (tile_id, ti, tj, tk).
	 */
	inline void global_coords_to_tile(
		uint32_t gi,
		uint32_t gj,
		uint32_t gk,
		uint32_t& tile_type,
		uint32_t& tile_id,
		uint32_t& ti,
		uint32_t& tj,
		uint32_t& tk
	) const
	{
		uint32_t tile_id_i = gi / m_tile_i_size;
		uint32_t tile_id_j = gj / m_tile_j_size;
		uint32_t tile_id_k = gk / m_tile_k_size;
		tile_id = tile_id_i * m_tile_j_num * m_tile_k_num + tile_id_j * m_tile_k_num + tile_id_k;
		tile_type = tile_id_to_tile_type(tile_id);

		ti = gi % m_tile_i_size;
		tj = gj % m_tile_j_size;
		tk = gk % m_tile_k_size;
	}

	/*
	 * Convert a cell from the tile coordinates from the perspective
	 * within a 32x32x32 tile (tile_id, ti, tj, tk) to its global
	 * coordinates from the perspective of the simulation box (gi, gj, gk).
	 */
	inline void tile_coords_to_global(
		uint32_t tile_id,
		uint32_t ti,
		uint32_t tj,
		uint32_t tk,
		uint32_t& gi,
		uint32_t& gj,
		uint32_t& gk
	) const
	{
		uint32_t tile_id_i = tile_id / (m_tile_j_num * m_tile_k_num);
		uint32_t tile_id_j = (tile_id - (tile_id_i * m_tile_j_num * m_tile_k_num)) / m_tile_k_num;
		uint32_t tile_id_k = tile_id - (tile_id_i * m_tile_j_num * m_tile_k_num + tile_id_j * m_tile_k_num);

		gi = tile_id_i * m_tile_i_size + ti;
		gj = tile_id_j * m_tile_j_size + tj;
		gk = tile_id_k * m_tile_k_size + tk;
	}

	/*
	 * Convert a cell from the tile coordinates from the perspective
	 * within a 32x32x32 tile (ti, tj, tk) to its subtile coordinates
	 * from the perspective of a 2x2x2 subtile (subtile_id, si, sj, sk).
	 */
	inline void regular_tile_coords_to_subtile(
		uint32_t ti,
		uint32_t tj,
		uint32_t tk,
		uint32_t& subtile_id,
		uint32_t& si,
		uint32_t& sj,
		uint32_t& sk
	) const
	{
		uint32_t subtile_id_i = ti / SUBTILE_I_SIZE;
		uint32_t subtile_id_j = tj / SUBTILE_J_SIZE;
		uint32_t subtile_id_k = tk / SUBTILE_K_SIZE;
		subtile_id = subtile_id_i * m_subtile_j_num * m_subtile_k_num + subtile_id_j * m_subtile_k_num + subtile_id_k;

		si = ti % SUBTILE_I_SIZE;
		sj = tj % SUBTILE_J_SIZE;
		sk = tk % SUBTILE_K_SIZE;
	}

	/*
	 * Convert a cell from the subtile coordinates from the perspective
	 * within a 2x2x2 subtile (subtile_id, si, sj, sk) to its tile
	 * coordinates from the perspective of a 32x32x32 tile (ti, tj, tk).
	 */
	inline void subtile_coords_to_regular_tile(
		uint32_t subtile_id,
		uint32_t si,
		uint32_t sj,
		uint32_t sk,
		uint32_t& ti,
		uint32_t& tj,
		uint32_t& tk
	) const
	{
		uint32_t subtile_id_i = subtile_id / (m_subtile_j_num * m_subtile_k_num);
		uint32_t subtile_id_j = (subtile_id - (subtile_id_i * m_subtile_j_num * m_subtile_k_num)) / m_subtile_k_num;
		uint32_t subtile_id_k = subtile_id - (subtile_id_i * m_subtile_j_num * m_subtile_k_num + subtile_id_j * m_subtile_k_num);

		ti = subtile_id_i * SUBTILE_I_SIZE + si;
		tj = subtile_id_j * SUBTILE_J_SIZE + sj;
		tk = subtile_id_k * SUBTILE_K_SIZE + sk;
	}

	/*
	 * Convert a cell from the tile coordinates from the perspective
	 * within a sparse tile (ti, tj, tk) to its subtile coordinates
	 * from the perspective of a 2x2x2 subtile (subtile_id, si, sj, sk).
	 */
	inline void sparse_tile_coords_to_subtile(
		uint32_t tile_type,
		uint32_t ti,
		uint32_t tj,
		uint32_t tk,
		uint32_t& subtile_id,
		uint32_t& si,
		uint32_t& sj,
		uint32_t& sk
	) const
	{
		uint32_t sparse_subtile_i_num = (tile_type & TILE_SPARSE_SUBTILE_I)
						? m_sparse_subtile_i_num
						: m_subtile_i_num;

		uint32_t sparse_subtile_j_num = (tile_type & TILE_SPARSE_SUBTILE_J)
						? m_sparse_subtile_j_num
						: m_subtile_j_num;

		uint32_t sparse_subtile_k_num = (tile_type & TILE_SPARSE_SUBTILE_K)
						? m_sparse_subtile_k_num
						: m_subtile_k_num;

		uint32_t subtile_id_i = ti / SUBTILE_I_SIZE;
		uint32_t subtile_id_j = tj / SUBTILE_J_SIZE;
		uint32_t subtile_id_k = tk / SUBTILE_K_SIZE;

		subtile_id = subtile_id_i * sparse_subtile_j_num * sparse_subtile_k_num;
		subtile_id += subtile_id_j * sparse_subtile_k_num;
		subtile_id += subtile_id_k;

		si = ti % SUBTILE_I_SIZE;
		sj = tj % SUBTILE_J_SIZE;
		sk = tk % SUBTILE_K_SIZE;
	}

	/*
	 * Convert a cell from the subtile coordinates from the perspective
	 * within a 2x2x2 subtile (subtile_id, si, sj, sk) to its tile
	 * coordinates from the perspective of a sparse tile (ti, tj, tk).
	 */
	inline void subtile_coords_to_sparse_tile(
		uint32_t subtile_id,
		uint32_t si,
		uint32_t sj,
		uint32_t sk,
		uint32_t tile_type,
		uint32_t& ti,
		uint32_t& tj,
		uint32_t& tk
	) const
	{
		uint32_t sparse_subtile_i_num = (tile_type & TILE_SPARSE_SUBTILE_I)
						? m_sparse_subtile_i_num
						: m_subtile_i_num;

		uint32_t sparse_subtile_j_num = (tile_type & TILE_SPARSE_SUBTILE_J)
						? m_sparse_subtile_j_num
						: m_subtile_j_num;

		uint32_t sparse_subtile_k_num = (tile_type & TILE_SPARSE_SUBTILE_K)
						? m_sparse_subtile_k_num
						: m_subtile_k_num;

		uint32_t subtile_id_i = subtile_id / (sparse_subtile_j_num * sparse_subtile_k_num);
		uint32_t subtile_id_j = (subtile_id - (subtile_id_i * sparse_subtile_j_num * sparse_subtile_k_num)) / sparse_subtile_k_num;
		uint32_t subtile_id_k = subtile_id - (subtile_id_i * sparse_subtile_j_num * sparse_subtile_k_num + subtile_id_j * sparse_subtile_k_num);

		ti = subtile_id_i * SUBTILE_I_SIZE + si;
		tj = subtile_id_j * SUBTILE_J_SIZE + sj;
		tk = subtile_id_k * SUBTILE_K_SIZE + sk;
	}

	inline void tile_coords_to_subtile(
		uint32_t tile_type,
		uint32_t ti,
		uint32_t tj,
		uint32_t tk,
		uint32_t& subtile_id,
		uint32_t& si,
		uint32_t& sj,
		uint32_t& sk
	) const
	{
		if (tile_type == TILE_REGULAR_SUBTILE) {
			return regular_tile_coords_to_subtile(ti, tj, tk, subtile_id, si, sj, sk);
		}
		else {
			return sparse_tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_id, si, sj, sk);
		}
	}

	inline void subtile_coords_to_tile(
		uint32_t subtile_id,
		uint32_t si,
		uint32_t sj,
		uint32_t sk,
		uint32_t tile_type,
		uint32_t& ti,
		uint32_t& tj,
		uint32_t& tk
	) const
	{
		if (tile_type == TILE_REGULAR_SUBTILE) {
			subtile_coords_to_regular_tile(subtile_id, si, sj, sk, ti, tj, tk);
		}
		else {
			subtile_coords_to_sparse_tile(subtile_id, si, sj, sk, tile_type, ti, tj, tk);
		}
	}

	/*
	 * Convert a cell from the subtile coordinates from the perspective
	 * within a 2x2x2 subtile (si, sj, sk) to its linear address from
	 * the perspective of a 1D array.
	 */
	inline void subtile_coords_to_linear(
		uint32_t n,
		uint32_t si,
		uint32_t sj,
		uint32_t sk,
		uint32_t& idx
	) const
	{
		const uint32_t n_size = 3;

		const uint32_t i_stride = SUBTILE_J_SIZE * SUBTILE_K_SIZE * n_size;
		const uint32_t j_stride = SUBTILE_K_SIZE * n_size;
		const uint32_t k_stride = n_size;

		idx = si * i_stride + sj * j_stride + sk * k_stride + n;
	}

	inline uint32_t tile_id_to_tile_type(uint32_t tile_id) const
	{
		uint32_t gi, gj, gk;
		tile_coords_to_global(
			tile_id, m_tile_i_size - 1, m_tile_j_size - 1, m_tile_k_size - 1,
			gi, gj, gk
		);

		uint32_t tile_type = TILE_REGULAR_SUBTILE;
		static_assert(TILE_REGULAR_SUBTILE == 0);

		tile_type |= (gi >= m_grid_unround_i_size) ? TILE_SPARSE_SUBTILE_I : 0;
		tile_type |= (gj >= m_grid_unround_j_size) ? TILE_SPARSE_SUBTILE_J : 0;
		tile_type |= (gk >= m_grid_unround_k_size) ? TILE_SPARSE_SUBTILE_K : 0;

		return tile_type;
	}
};

/*
 * GlobalArray is a wrapper around a Kokkos's multi-dimensional view.
 * It's used to store both fields and operators of the entire simulation
 * box. Preferable, it should be accessed at the granularity of an
 * entire tile for memory spatial locality - which is usually the case
 * during a simulation (note that Kokkos also automatically performs
 * initialization with proper NUMA first-touch, at the outermost dimension).
 */
template <typename T>
struct KokkosTile;	/* incomplete type, to be declared later */

template <typename T>
struct KokkosGlobalArray
{
	KokkosGlobalArray(
		const std::string& name,
		KokkosGrid& grid
	) :
		m_grid(grid)
	{

		m_view = new Kokkos::View<
			KokkosSubtile<T>**,
			Kokkos::DefaultExecutionSpace::memory_space,
			Kokkos::MemoryTraits<Kokkos::Restrict>
		>(name, grid.m_tile_num, grid.m_subtile_num);
	}

	inline const KokkosTile<T> get_tile(uint32_t tile_id) const
	{
		return KokkosTile<T>(m_grid, *this, tile_id);
	}

#if 0
	inline const KokkosTile<T> get_tile(
		uint32_t gi,
		uint32_t gj,
		uint32_t gk
	) const
	{
		uint32_t tile_id, ti, tj, tk;

		grid.global_coords_to_tile(gi, gj, gk, tile_id, ti, tj, tk);
		return KokkosTile<T>(m_grid, *this, tile_id, TILE_SPARSE_SUBTILE);
	}


	inline T& operator() (
		uint32_t tile_id,
		uint32_t n,
		uint32_t ti,
		uint32_t tj,
		uint32_t tk
	) const
	{
		uint32_t subtile_id, si, sj, sk;
		uint32_t idx;

                m_grid.tile_coords_to_subtile(ti, tj, tk, subtile_id, si, sj, sk);
                m_grid.subtile_coords_to_linear(n, si, sj, sk, idx);
		
		return (*m_view)(tile_id, subtile_id).array[idx]; 
	}
#endif

	inline T& operator() (
		uint32_t n,
		uint32_t gi,
		uint32_t gj,
		uint32_t gk
	)
	{
		uint32_t tile_type, tile_id, ti, tj, tk;
		uint32_t subtile_id, si, sj, sk;
		uint32_t idx;
		m_grid.global_coords_to_tile(gi, gj, gk, tile_type, tile_id, ti, tj, tk);
		m_grid.tile_coords_to_subtile(tile_type, ti, tj, tk, subtile_id, si, sj, sk);
		m_grid.subtile_coords_to_linear(n, si, sj, sk, idx);

		return get_tile(tile_id).get_subtile(subtile_id).array[idx];
	}

	KokkosGrid& m_grid;
	Kokkos::View<
		KokkosSubtile<T>**,
		Kokkos::DefaultExecutionSpace::memory_space,
		Kokkos::MemoryTraits<Kokkos::Restrict>
	>* __restrict__ m_view;
};

/*
 * Tile is a wrapper of the same underlying memory that belongs to a
 * GlobalArray, but provides an accessor for addressing it at the
 * granularity of a tile.
 */
template <typename T>
struct KokkosTile
{
	KokkosTile(
		KokkosGrid& grid,
		const KokkosGlobalArray<T>& array,
		uint32_t tile_id
	) :
		m_grid(grid),
		m_array(array),
		m_tile_id(tile_id)
	{};

	inline KokkosSubtile<T>& get_subtile(uint32_t subtile_id) const
	{
		return (*m_array.m_view)(m_tile_id, subtile_id);
	}

	KokkosGrid& m_grid;
	const KokkosGlobalArray<T>& m_array;
	uint32_t m_tile_id;
};


/*
 * LocalTile is a wrapper around a Kokkos's per-team memory, which is a
 * temporary buffer. On the GPU, it corresponds to a GPU's on-chip memory
 * shared by a workgroup (not implemented in this engine). On the CPU,
 * it's just regular thread-local memory, but its size is often set equal
 * to a CPU's cache size to express cache locality in the code explicitly.
 *
 * In the engine, loading and saving the current tile from/to GlobalArray
 * is the beginning and end of each loop.
 */
template <typename T>
struct KokkosLocalTile
{
	KokkosLocalTile(
		const KokkosGrid& grid,
		const Kokkos::TeamPolicy<>::member_type& team_member
	) :
		m_grid(grid)
	{
		auto& shmem = team_member.team_shmem();
		
		size_t shmem_bytes = grid.m_subtile_num * sizeof(KokkosSubtile<T>);
		m_shmem_ptr = (KokkosSubtile<T>*) shmem.get_shmem(shmem_bytes);

		if (!m_shmem_ptr) {
			fprintf(stderr, "KokkosLocalTile: get_shmem() failed!\n");
			abort();
		}
	}

	const KokkosGrid& m_grid;
	KokkosSubtile<T>* __restrict__ m_shmem_ptr;
	uint32_t m_type;
	uint32_t subtile_num;
	uint32_t tile_i_size, tile_j_size, tile_k_size;

#if 0
	inline T& operator() (
		uint32_t n,
		uint32_t ti,
		uint32_t tj,
		uint32_t tk
	)
	{
		uint32_t subtile_id, si, sj, sk;
		uint32_t idx;

		m_grid.tile_coords_to_subtile(m_type, ti, tj, tk, subtile_id, si, sj, sk);
                m_grid.subtile_coords_to_linear(n, si, sj, sk, idx);
		
		return m_shmem_ptr[subtile_id].array[idx];
	}
#endif

	inline KokkosSubtile<T>& get_subtile(uint32_t subtile_id) const
	{
		return m_shmem_ptr[subtile_id];
	}

	void load_from(
		const uint32_t tile_id,
		const KokkosTile<T>& tile
	)
	{
		m_type = m_grid.tile_id_to_tile_type(tile_id);

		tile_i_size = (m_type & TILE_SPARSE_SUBTILE_I)
				? m_grid.m_sparse_tile_i_size
				: m_grid.m_tile_i_size;

		tile_j_size = (m_type & TILE_SPARSE_SUBTILE_J)
				? m_grid.m_sparse_tile_j_size
				: m_grid.m_tile_j_size;

		tile_k_size = (m_type & TILE_SPARSE_SUBTILE_K)
				? m_grid.m_sparse_tile_k_size
				: m_grid.m_tile_k_size;

		if (m_type == TILE_REGULAR_SUBTILE) {
			subtile_num = m_grid.m_subtile_num;
		}
		else {
			uint32_t sparse_subtile_i_num = (m_type & TILE_SPARSE_SUBTILE_I)
							? m_grid.m_sparse_subtile_i_num
							: m_grid.m_subtile_i_num;

			uint32_t sparse_subtile_j_num = (m_type & TILE_SPARSE_SUBTILE_J)
							? m_grid.m_sparse_subtile_j_num
							: m_grid.m_subtile_j_num;

			uint32_t sparse_subtile_k_num = (m_type & TILE_SPARSE_SUBTILE_K)
							? m_grid.m_sparse_subtile_k_num
							: m_grid.m_subtile_k_num;

			subtile_num = sparse_subtile_i_num * sparse_subtile_j_num * sparse_subtile_k_num;
		}

		for (uint32_t subtile_id = 0; subtile_id < subtile_num; subtile_id++) {
			const KokkosSubtile<T>& subtile = tile.get_subtile(subtile_id);
			for (uint32_t i = 0; i < SUBTILE_SIZE * 3; i++) {
				m_shmem_ptr[subtile_id].array[i] = subtile.array[i];
				//fprintf(stderr, "copy %zu to %zu\n", &subtile.array[i], &(m_shmem_ptr[subtile_id].array[i]));
			}
		}
	}

	void save_to(
		const uint32_t tile_id,
		const KokkosTile<T>& tile
	) const
	{
		for (uint32_t subtile_id = 0; subtile_id < subtile_num; subtile_id++) {
			KokkosSubtile<T>& subtile = tile.get_subtile(subtile_id);
			memcpy(subtile.array, m_shmem_ptr[subtile_id].array, SUBTILE_SIZE * 3 * sizeof(T));
			//for (uint32_t i = 0; i < SUBTILE_SIZE * 3; i++) {
			//	subtile.array[i] = m_shmem_ptr[subtile_id].array[i];
			//}
		}
	}
};
