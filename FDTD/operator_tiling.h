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

#ifndef OPERATOR_TILING_H
#define OPERATOR_TILING_H

#include "operator_multithread.h"

#include <boost/thread.hpp>

#define OPERATOR_TILING_BASE Operator_Multithread

class Operator_Tiling : public OPERATOR_TILING_BASE
{
	friend class Engine_Multithread;
	friend class Engine_Tiling;
	friend class Operator_Thread;
public:
	//! Create a new operator
	static Operator_Tiling* New(unsigned int numThreads = 0);
	virtual Engine* CreateEngine();
};

#endif // OPERATOR_TILING_H
