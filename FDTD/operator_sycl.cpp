/*
*	Copyright (C) 2010 Sebastian Held (Sebastian.Held@gmx.de)
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

#include "engine_sycl.h"
#include "operator_sycl.h"
#include "tools/array_ops.h"

Operator_sycl* Operator_sycl::New()
{
	cout << "Create GPU-accelerated FDTD operator (SYCL)" << endl;
	Operator_sycl* op = new Operator_sycl();
	op->Init();
	return op;
}

Operator_sycl::Operator_sycl() : Operator()
{
	f4_vv_ptr = 0;
	f4_vi_ptr = 0;
	f4_iv_ptr = 0;
	f4_ii_ptr = 0;
}

Operator_sycl::~Operator_sycl()
{
	Delete();
}

Engine* Operator_sycl::CreateEngine()
{
	//! create a special sse-engine
	m_Engine = Engine_sycl::New(this);
	return m_Engine;
}

void Operator_sycl::Init()
{
	Operator::Init();
	f4_vv_ptr = 0;
	f4_vi_ptr = 0;
	f4_iv_ptr = 0;
	f4_ii_ptr = 0;
}

void Operator_sycl::Delete()
{
	Delete_SYCL_N_3DArray(f4_vv_ptr,numLines);
	Delete_SYCL_N_3DArray(f4_vi_ptr,numLines);
	Delete_SYCL_N_3DArray(f4_iv_ptr,numLines);
	Delete_SYCL_N_3DArray(f4_ii_ptr,numLines);
	f4_vv_ptr = 0;
	f4_vi_ptr = 0;
	f4_iv_ptr = 0;
	f4_ii_ptr = 0;
}

void Operator_sycl::Reset()
{
	Delete();
	Operator::Reset();
}


void Operator_sycl::InitOperator()
{
	Delete_SYCL_N_3DArray(f4_vv_ptr,numLines);
	Delete_SYCL_N_3DArray(f4_vi_ptr,numLines);
	Delete_SYCL_N_3DArray(f4_iv_ptr,numLines);
	Delete_SYCL_N_3DArray(f4_ii_ptr,numLines);
	f4_vv_ptr = Create_SYCL_N_3DArray<sycl::float4>(m_sycl_queue, numLines);
	f4_vi_ptr = Create_SYCL_N_3DArray<sycl::float4>(m_sycl_queue, numLines);
	f4_iv_ptr = Create_SYCL_N_3DArray<sycl::float4>(m_sycl_queue, numLines);
	f4_ii_ptr = Create_SYCL_N_3DArray<sycl::float4>(m_sycl_queue, numLines);

	numVectors =  ceil((double)numLines[2]/4.0);
}
