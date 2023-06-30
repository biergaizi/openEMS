/*
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

#ifndef ENGINE_MULTITHREAD_H
#define ENGINE_MULTITHREAD_H

#include "operator_multithread.h"
#include "engine_sse_compressed.h"
#include "tools/tiling.h"

#include <boost/thread.hpp>
#include <boost/fusion/include/list.hpp>
#include <boost/fusion/container/list/list_fwd.hpp>
#include <boost/fusion/include/list_fwd.hpp>

#include "tools/useful.h"
#ifndef __GNUC__
#include <Winsock2.h> // for struct timeval
#else
#include <sys/time.h>
#endif


#ifdef MPI_SUPPORT
	#define ENGINE_MULTITHREAD_BASE Engine_MPI
	#include "engine_mpi.h"
#else
	#define ENGINE_MULTITHREAD_BASE Engine_SSE_Compressed
#endif


class Engine_Multithread;

namespace NS_Engine_Multithread
{

class DBG   // debug
{
public:
	DBG() {}
	~DBG() { std::cout << os.str();}
	std::ostringstream& cout() {return os;}
protected:
	std::ostringstream os;
};

class Timer   //debug
{
public:
	Timer() {gettimeofday(&t1,NULL);}
	double elapsed() {gettimeofday(&t2,NULL); return (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)*1e-6;}
protected:
	timeval t1,t2;
};

class thread
{
public:
	thread( Engine_Multithread* ptr, std::vector<Tiles3D> tiles, unsigned int blkTimesteps, std::vector<Tiles3D> fallbackTiles, unsigned int threadID );
	void operator()();
protected:
	unsigned int m_start, m_stop, m_stop_h, m_threadID;
	unsigned int m_blkTimesteps;
	std::vector<Tiles3D> m_tiles, m_fallbackTiles;
	Engine_Multithread *m_enginePtr;
	void iterateTimesteps(int timestep, std::vector<Range3D>&);
	void iterateUnskewedSingleTimestep(int timestep, std::vector<Range3D>&);
};
} // namespace


class Engine_Multithread : public ENGINE_MULTITHREAD_BASE
{
	friend class NS_Engine_Multithread::thread;
	friend class Engine_CylinderMultiGrid;
public:
	static Engine_Multithread* New(const Operator_Multithread* op, unsigned int numThreads = 0);
	virtual ~Engine_Multithread();

	virtual void setNumThreads( unsigned int numThreads );
	virtual void Init();
	virtual void Reset();
	virtual void NextInterval(float curr_speed);

	//! Iterate \a iterTS number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual void DoPreVoltageUpdates(int threadID, int timestep, int start[3], int end[3]);
	virtual void DoPostVoltageUpdates(int threadID, int timestep, int start[3], int end[3]);
	virtual void Apply2Voltages(int threadID, int timestep, int start[3], int end[3]);

	virtual void DoPreCurrentUpdates(int threadID, int timestep, int start[3], int end[3]);
	virtual void DoPostCurrentUpdates(int threadID, int timestep, int start[3], int end[3]);
	virtual void Apply2Current(int threadID, int timestep, int start[3], int end[3]);

	virtual void InitializeTiling(std::vector<Range3D> tiles);

protected:
	Engine_Multithread(const Operator_Multithread* op);
	void changeNumThreads(unsigned int numThreads);
	const Operator_Multithread* m_Op_MT;
	boost::thread_group *m_thread_group;
	boost::barrier *m_startBarrier, *m_stopBarrier;
	boost::barrier *m_IterateBarrier;
	volatile unsigned int m_iterTS;
	unsigned int m_numThreads; //!< number of worker threads
	unsigned int m_max_numThreads; //!< max. number of worker threads
	volatile bool m_stopThreads;
	bool m_opt_speed;
	float m_last_speed;

#ifdef MPI_SUPPORT
	/*! Workaround needed for subgridding scheme... (see Engine_CylinderMultiGrid)
	 Some engines may need an additional barrier for synchronizing MPI communication.
	 This engine will not initialize or cleanup this barrier, but check for it and wait before executing any MPI sync.
	 Make sure to cleanup (delete) this barriere before Engine_Multithread::Reset() is called.
	 */
	boost::barrier *m_MPI_Barrier;
#endif

#ifdef ENABLE_DEBUG_TIME
	std::map<boost::thread::id, std::vector<double> > m_timer_list;
#endif
};

#endif // ENGINE_MULTITHREAD_H
