/*
*	Copyright (C) 2023 Yifeng Li (tomli@tomli.me)
*	Copyright (C) 2010 Sebastian Held (sebastian.held@gmx.de)
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

#include "engine_tiling.h"
#include "extensions/engine_extension.h"
#include "tools/array_ops.h"

#include <iomanip>

#ifndef SSE_CORRECT_DENORMALS
#include <xmmintrin.h>
#endif

#include "tools/tiling.h"

Engine_Tiling* Engine_Tiling::New(const Operator_Multithread* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (compressed SSE + multi-threading + spatial/temporal tiling)" << endl;
	Engine_Tiling* e = new Engine_Tiling(op);
	e->setNumThreads( numThreads );
	e->Init();
	return e;
}

Engine_Tiling::Engine_Tiling(const Operator_Multithread* op) : ENGINE_TILING_BASE(op)
{
	m_Op_MT = op;
	m_type = SSE;
	m_IterateBarrier = 0;
	m_startBarrier = 0;
	m_stopBarrier = 0;
	m_thread_group = 0;
	m_max_numThreads = boost::thread::hardware_concurrency();
	m_numThreads = 0;
	m_last_speed = 0;
	m_opt_speed = false;
	m_stopThreads = true;
}

Engine_Tiling::~Engine_Tiling()
{
#ifdef MPI_SUPPORT
	std::cerr << "Tiling engine does not support MPI!" << std::endl;
	std::exit(1);
#endif

	Reset();
}

void Engine_Tiling::setNumThreads( unsigned int numThreads )
{
	m_numThreads = numThreads;
}

void Engine_Tiling::Init()
{
	m_stopThreads = true;
	m_opt_speed = false;
	ENGINE_TILING_BASE::Init();

	// initialize threads
	m_stopThreads = false;
	if (m_numThreads == 0)
	{
		m_opt_speed = true;
		m_numThreads = 1;
	}
	//else if (m_numThreads > m_max_numThreads)
	//	m_numThreads = m_max_numThreads;

	this->changeNumThreads(m_numThreads);
}

void Engine_Tiling::Reset()
{
	if (m_thread_group!=0) // prevent multiple invocations
	{
		ClearExtensions(); //prevent extensions from interfering with thread reset...

		// stop the threads
		//NS_Engine_Tiling::DBG().cout() << "stopping all threads" << endl;
		m_thread_group->interrupt_all();
		m_thread_group->join_all(); // wait for termination
		delete m_IterateBarrier;
		m_IterateBarrier = 0;
		delete m_startBarrier;
		m_startBarrier = 0;
		delete m_stopBarrier;
		m_stopBarrier = 0;
		delete m_thread_group;
		m_thread_group = 0;
	}

	ENGINE_TILING_BASE::Reset();
}

void Engine_Tiling::changeNumThreads(unsigned int numThreads)
{
	if (m_thread_group!=0)
	{
		m_thread_group->interrupt_all();
		m_thread_group->join_all(); // wait for termination
		delete m_thread_group;
		m_thread_group = 0;
		//m_stopThreads = false;
	}

	m_numThreads = numThreads;

	if (g_settings.GetVerboseLevel()>0)
		cout << "Tilinged engine using " << m_numThreads << " threads. Utilization: (";

	vector<unsigned int> m_Start_Lines;
	vector<unsigned int> m_Stop_Lines;
	m_Op_MT->CalcStartStopLines( m_numThreads, m_Start_Lines, m_Stop_Lines );

	if (m_IterateBarrier!=0)
		delete m_IterateBarrier;
	// make sure all threads are waiting
	m_IterateBarrier = new boost::barrier(m_numThreads); // numThread workers

	if (m_startBarrier!=0)
		delete m_startBarrier;
	m_startBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller

	if (m_stopBarrier!=0)
		delete m_stopBarrier;
	m_stopBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller

	// TODO: all the casting between "unsigned int" and "int" should be eliminated by
	// adjusting tiling.cpp to use "unsigned int", which is currently not done since
	// tiling.cpp uses "-1" for invalid values internally.
	int blkSizes[3] = {
		//m_Op_MT->GetNumberOfLines(0),
		//m_Op_MT->GetNumberOfLines(1),
		10,
		10,
		(int) m_Op_MT->GetNumberOfLines(2)
	};
	int blkHalfTimesteps = 5 * 2;

	fprintf(stderr, "calculating tiling for X, Y and Z axis for %d timesteps\n", blkHalfTimesteps / 2);
	Tiles tilesX = computeDiamondTiles1D(numLines[0], blkSizes[0], blkHalfTimesteps);
	Tiles tilesY = computeDiamondTiles1D(numLines[1], blkSizes[1], blkHalfTimesteps);
	Tiles tilesZ = computeDiamondTiles1D(numLines[2], blkSizes[2], blkHalfTimesteps);

	fprintf(stderr, "combining tilings from X, Y and Z axis\n");
        auto tilesPerStagePerThread = combineTilesTo3D(
		tilesX, tilesY, tilesZ,
		blkHalfTimesteps,
		m_numThreads
	);

	fprintf(stderr, "flattening tiling\n");
	std::vector<Range3D> tilesForAllThreads;
	for (auto& tilesPerStage : tilesPerStagePerThread)
	{
		for (auto& tiles : tilesPerStage)
		{
			for (auto& tile : tiles)
			{
				tilesForAllThreads.push_back(tile);
			}
		}
	}

	// This is rather ugly but we do need a fallback without temporal tiling,
	// but at least we can still use spatial tiling.
	auto fallbackTilesPerStagePerThread = computeRectangularTiles3D((int*) numLines, blkSizes, m_numThreads);
	for (auto& fallbackTilesPerStage : fallbackTilesPerStagePerThread)
	{
		for (auto& tiles : fallbackTilesPerStage)
		{
			for (auto& tile : tiles)
			{
				tilesForAllThreads.push_back(tile);
			}
		}
	}

	m_thread_group = new boost::thread_group();
	for (unsigned int n=0; n<m_numThreads; n++)
	{
		boost::thread *t = new boost::thread(
			NS_Engine_Tiling::thread(
				this,
				tilesPerStagePerThread[n], blkHalfTimesteps / 2,
				fallbackTilesPerStagePerThread[n],
				n
		));
		m_thread_group->add_thread( t );
	}

	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->SetNumberOfThreads(m_numThreads);

	fprintf(stderr, "hashing tiling\n");
	InitializeTiling(tilesForAllThreads);
}

bool Engine_Tiling::IterateTS(unsigned int iterTS)
{
	m_iterTS = iterTS;
	//fprintf(stderr, "%d TS requested.", iterTS);

	m_startBarrier->wait(); // start the threads
	m_stopBarrier->wait(); // wait for the threads to finish <iterTS> time steps

	return true;
}

void Engine_Tiling::NextInterval(float curr_speed)
{
	ENGINE_TILING_BASE::NextInterval(curr_speed);
	if (!m_opt_speed) return;
	if (curr_speed<m_last_speed)
	{
		this->changeNumThreads(m_numThreads-1);
		cout << "Tilinged Engine: Best performance found using " << m_numThreads << " threads." << std::endl;
		m_opt_speed = false;
	}
	else if (m_numThreads<m_max_numThreads)
	{
		m_last_speed = curr_speed;
		this->changeNumThreads(m_numThreads+1);
	}
}

void Engine_Tiling::DoPreVoltageUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in reverse order -> highest priority gets access to the voltages last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
	{
		m_Eng_exts.at(n)->DoPreVoltageUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}

}

void Engine_Tiling::DoPostVoltageUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->DoPostVoltageUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::Apply2Voltages(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->Apply2Voltages(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::DoPreCurrentUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in reverse order -> highest priority gets access to the currents last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
	{
		m_Eng_exts.at(n)->DoPreCurrentUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::DoPostCurrentUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->DoPostCurrentUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::Apply2Current(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->Apply2Current(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::InitializeTiling(std::vector<Range3D> tiles)
{
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		std::cerr << "Using ";
		std::cerr << m_Eng_exts.at(n)->GetExtensionName();
		std::cerr << std::endl;

		m_Eng_exts.at(n)->InitializeTiling(tiles);
	}
}

//
// *************************************************************************************************************************
//
namespace NS_Engine_Tiling
{

thread::thread(
	Engine_Tiling* ptr,
	std::vector<Tiles3D> tiles,
	unsigned int blkTimesteps,
	std::vector<Tiles3D> fallbackTiles,
	unsigned int threadID
)
{
	m_enginePtr = ptr;
	m_tiles = tiles;
	m_threadID = threadID;
	m_fallbackTiles = fallbackTiles;
	m_blkTimesteps = blkTimesteps;

	fprintf(stderr, "Tiles - Thread: %d / Tiles: %ld\n", threadID, tiles.size());
}

void thread::operator()()
{
	// speed up the calculation of denormal floating point values (flush-to-zero)
#ifndef SSE_CORRECT_DENORMALS
	unsigned int oldMXCSR = _mm_getcsr(); //read the old MXCSR setting
	unsigned int newMXCSR = oldMXCSR | 0x8040; // set DAZ and FZ bits
	_mm_setcsr( newMXCSR ); //write the new MXCSR setting to the MXCSR
#endif

	//fprintf(stderr, "Tiles - X: %d / Y: %d / Z: %d\n", blk_x_max, blk_y_max, blk_z_max);

	while (!m_enginePtr->m_stopThreads)
	{
		// wait for start
		m_enginePtr->m_startBarrier->wait();

		if (m_enginePtr->m_stopThreads)
			return;

		// how many timesteps are we calculating at the same time in
		// diamond tiling?
		unsigned int batchTimesteps = m_enginePtr->m_iterTS / m_blkTimesteps;
		unsigned int leftoverTimesteps = m_enginePtr->m_iterTS % m_blkTimesteps;

		int currentTimestep = m_enginePtr->numTS;

		for (unsigned int iter=0; iter < batchTimesteps; ++iter)
		{
			for (unsigned int stage = 0; stage < m_tiles.size(); stage++)
			{
				iterateTimesteps(currentTimestep, m_tiles[stage]);
				m_enginePtr->m_IterateBarrier->wait();
			}
			//std::exit(0);

			currentTimestep += m_blkTimesteps;
		}

		if (m_threadID == 0 && leftoverTimesteps > 0) {
			//fprintf(stderr, "requested %d TS\n", m_enginePtr->m_iterTS);
			//fprintf(stderr, "this iteration has %d leftover timesteps...\n", leftoverTimesteps);
		}

		for (unsigned int iter = 0; iter < leftoverTimesteps; ++iter)
		{
			for (unsigned int stage = 0; stage < m_fallbackTiles.size(); stage++)
			{
				iterateUnskewedSingleTimestep(currentTimestep, m_fallbackTiles[stage]);
			}
			currentTimestep += 1;

		}

		if (m_threadID == 0) {
			 // only the first thread increments numTS
			m_enginePtr->numTS = currentTimestep;
		}
		m_enginePtr->m_stopBarrier->wait();
	}
}

void thread::iterateTimesteps(int timestep, std::vector<Range3D>& tiles)
{
	auto op = m_enginePtr->m_Op_MT;
	int baseTimestep = timestep;

	for (auto& tile : tiles)
	{
		int timestep = baseTimestep + tile.timestep;

		// pre voltage stuff...
		m_enginePtr->DoPreVoltageUpdates(timestep, tile.voltageStart, tile.voltageStop);

		m_enginePtr->UpdateVoltages(
			(unsigned int *) tile.voltageStart,
			(unsigned int *) tile.voltageStop
		);

		//post voltage stuff...
		m_enginePtr->DoPostVoltageUpdates(timestep, tile.voltageStart, tile.voltageStop);
		m_enginePtr->Apply2Voltages(timestep, tile.voltageStart, tile.voltageStop);

		//pre current stuff
		m_enginePtr->DoPreCurrentUpdates(timestep, tile.currentStart, tile.currentStop);

		unsigned int currentStopSkipLast[3] = {
			(unsigned int) tile.currentStop[0],
			(unsigned int) tile.currentStop[1],
			(unsigned int) tile.currentStop[2]
		};

		for (int n = 0; n < 3; n++)
		{
			if (currentStopSkipLast[n] > op->GetNumberOfLines(n) - 2)
			{
				currentStopSkipLast[n] = op->GetNumberOfLines(n) - 2;
			}
		}

		//current updates
		m_enginePtr->UpdateCurrents((unsigned int *) tile.currentStart, currentStopSkipLast);

		//post current stuff
		m_enginePtr->DoPostCurrentUpdates(timestep, tile.currentStart, tile.currentStop);
		m_enginePtr->Apply2Current(timestep, tile.currentStart, tile.currentStop);
	}
}

void thread::iterateUnskewedSingleTimestep(int timestep, std::vector<Range3D>& tiles)
{
	auto op = m_enginePtr->m_Op_MT;

	for (auto& tile : tiles)
	{
		// pre voltage stuff...
		m_enginePtr->DoPreVoltageUpdates(timestep, tile.voltageStart, tile.voltageStop);

		m_enginePtr->UpdateVoltages(
			(unsigned int *) tile.voltageStart,
			(unsigned int *) tile.voltageStop
		);

		//post voltage stuff...
		m_enginePtr->DoPostVoltageUpdates(timestep, tile.voltageStart, tile.voltageStop);
		m_enginePtr->Apply2Voltages(timestep, tile.voltageStart, tile.voltageStop);
	}

	m_enginePtr->m_IterateBarrier->wait();

	for (auto& tile : tiles) {
		//pre current stuff
		m_enginePtr->DoPreCurrentUpdates(timestep, tile.currentStart, tile.currentStop);

		unsigned int currentStopSkipLast[3] = {
			(unsigned int) tile.currentStop[0],
			(unsigned int) tile.currentStop[1],
			(unsigned int) tile.currentStop[2]
		};

		for (int n = 0; n < 3; n++)
		{
			if (currentStopSkipLast[n] > op->GetNumberOfLines(n) - 2)
			{
				currentStopSkipLast[n] = op->GetNumberOfLines(n) - 2;
			}
		}

		//current updates
		m_enginePtr->UpdateCurrents((unsigned int *) tile.currentStart, currentStopSkipLast);

		//post current stuff
		m_enginePtr->DoPostCurrentUpdates(timestep, tile.currentStart, tile.currentStop);
		m_enginePtr->Apply2Current(timestep, tile.currentStart, tile.currentStop);
	}

	m_enginePtr->m_IterateBarrier->wait();
}

} // namespace
