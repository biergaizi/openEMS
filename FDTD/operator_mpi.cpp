/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY{} without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "operator_mpi.h"
#include "operator_sse_compressed.h"
#include "engine_sse_compressed.h"
#include "engine_mpi.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "mpi.h"

Operator_MPI* Operator_MPI::New()
{
	cout << "Create FDTD operator (compressed SSE + MPI)" << endl;
	Operator_MPI* op = new Operator_MPI();
	op->Init();
	return op;
}

Operator_MPI::Operator_MPI() : Operator_SSE_Compressed()
{
	m_NumProc = MPI::COMM_WORLD.Get_size();

	//enabled only if more than one process is active
	m_MPI_Enabled = m_NumProc>0;
}

Operator_MPI::~Operator_MPI()
{
	Delete();
}

bool Operator_MPI::SetGeometryCSX(ContinuousStructure* geo)
{
	//manipulate geometry for this part...

	if (m_MPI_Enabled)
	{
		CSRectGrid* grid = geo->GetGrid();
		int nz = grid->GetQtyLines(2);
		std::vector<unsigned int> jobs = AssignJobs2Threads(nz, m_NumProc);
		double z_lines[jobs.at(m_MyID)+1];

		if (m_MyID==0)
		{
			for (unsigned int n=0;n<jobs.at(0);++n)
				z_lines[n] = grid->GetLine(2,n);
			grid->ClearLines(2);
			grid->AddDiscLines(2,jobs.at(0),z_lines);

		}
		else
		{
			unsigned int z_start=0;
			for (int n=0;n<m_MyID;++n)
				z_start+=jobs.at(n);
			for (unsigned int n=0;n<=jobs.at(m_MyID);++n)
				z_lines[n] = grid->GetLine(2,z_start+n-1);
			grid->ClearLines(2);
			grid->AddDiscLines(2,jobs.at(m_MyID)+1,z_lines);
		}

		//lower neighbor is ID-1
		if (m_MyID>0)
			m_NeighborDown[2]=m_MyID-1;
		//upper neighbor is ID+1
		if (m_MyID<m_NumProc-1)
			m_NeighborUp[2]=m_MyID+1;
	}
	else
		cerr << "Operator_MPI::SetGeometryCSX: Warning: Number of MPI processes is 1, skipping MPI engine... " << endl;

	return Operator_SSE_Compressed::SetGeometryCSX(geo);
}

double Operator_MPI::CalcTimestep()
{
	double ret = Operator::CalcTimestep();

	if (!m_MPI_Enabled)
		return ret;

	double local_dT = dT;
	//find the smallest time-step requestes by all processings
	MPI_Reduce(&local_dT, &dT, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	//send the smallest time-step to all
	MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	return ret;
}


void Operator_MPI::SetBoundaryCondition(int* BCs)
{
	if (!m_MPI_Enabled)
		return Operator_SSE_Compressed::SetBoundaryCondition(BCs);

	//set boundary conditions on MPI interfaces to PEC, ApplyElectricBC will handle proper interface handling...
	for (int n=0;n<3;++n)
	{
		if (m_NeighborUp[n]>=0)
			BCs[2*n+1] = 0;
		if (m_NeighborDown[n]>=0)
			BCs[2*n] = 0;
	}
	Operator_SSE_Compressed::SetBoundaryCondition(BCs);
}

void Operator_MPI::ApplyElectricBC(bool* dirs)
{
	if (!m_MPI_Enabled)
		return Operator_SSE_Compressed::ApplyElectricBC(dirs);

	for (int n=0;n<3;++n)
	{
		//do not delete operator at upper inteface
		if (m_NeighborUp[n]>=0)
			dirs[2*n+1] = false;
	}
	Operator_SSE_Compressed::ApplyElectricBC(dirs);
}

Engine* Operator_MPI::CreateEngine() const
{
	if (m_MPI_Enabled)
		return Engine_MPI::New(this);
	else
		return Engine_SSE_Compressed::New(this);
}

void Operator_MPI::Init()
{
	Operator_SSE_Compressed::Init();

	m_MyTag = 0;

	for (int i=0;i<3;++i)
	{
		m_NeighborUp[i]=-1;
		m_NeighborDown[i]=-1;
	}

	int  namelen;
	m_NumProc = MPI::COMM_WORLD.Get_size();
	m_MyID = MPI::COMM_WORLD.Get_rank();

	m_Processor_Name = new char[MPI_MAX_PROCESSOR_NAME];
	MPI::Get_processor_name(m_Processor_Name,namelen);

	if (m_MPI_Enabled)
		cerr << "Operator_MPI::Init(): Running on " << m_Processor_Name << endl;
}

void Operator_MPI::Delete()
{
	delete[] m_Processor_Name;
}

void Operator_MPI::Reset()
{
	Delete();
	Operator_SSE_Compressed::Reset();
}

string Operator_MPI::PrependRank(string name)
{
	stringstream out_name;
	if (m_MPI_Enabled)
		out_name << "ID" << m_MyID << "_" << name;
	else
		out_name << name;
	return out_name.str();
}

void Operator_MPI::DumpOperator2File(string filename)
{
	Operator_SSE_Compressed::DumpOperator2File(PrependRank(filename));
}

void Operator_MPI::DumpMaterial2File(string filename)
{
	Operator_SSE_Compressed::DumpMaterial2File(PrependRank(filename));
}

void Operator_MPI::DumpPEC2File( string filename )
{
	Operator_SSE_Compressed::DumpPEC2File(PrependRank(filename));
}
