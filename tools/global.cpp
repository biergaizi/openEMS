/*
*	Copyright (C) 2010 Sebastian Held <sebastian.held@gmx.de>
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

#include <cstring>
#include <iostream>
#include "global.h"

using namespace std;

// create global object
Global g_settings;

Global::Global()
{
	m_showProbeDiscretization = false;
	m_nativeFieldDumps = false;
	m_VerboseLevel = 0;
}

namespace po = boost::program_options;

po::options_description
Global::cmdArgs()
{
	po::options_description optdesc("Additional global arguments");
	optdesc.add_options()
		("showProbeDiscretization", "Show probe discretization information")
		("nativeFieldDumps",        "Dump all fields using the native field "
		                            "components")
		("verbose",                 po::value<unsigned int>()->default_value(0),
									"Verbose level, select debug level 1 to 3, "
		                            "also accept -v, -vv, -vvv");
	return optdesc;
}

//! \brief This function initializes the object
void Global::parseCommandLineArguments(int argc, const char* argv[])
{
	// Hack: boost::program_options doesn't support repeated "-vv"
	// and "-vvv" syntax, it also doesn't support omitting the
	// value of a parameter, both causes validation failure. It's
	// not worthwhile to write a custom validator just for exactly
	// a single special case. Just change argv[] to avoid them.
	// This must be the first parseCommandLineArguments() to call
	// in main().
	for (int i = 0; i < argc; i++)
	{
		const char* arg = argv[i];

		if (strcmp(arg, "--verbose") == 0)
		{
			argv[i] = "--verbose=1";
		}
		if (strcmp(arg, "-v") == 0)
		{
			argv[i] = "--verbose=1";
		}
		if (strcmp(arg, "-vv") == 0)
		{
			argv[i] = "--verbose=2";
		}
		else if (strcmp(arg, "-vvv") == 0)
		{
			argv[i] = "--verbose=3";
		}
	}

	po::options_description opts = cmdArgs();
	po::variables_map map;
	po::store(po::command_line_parser(argc, argv).options(opts).
		allow_unregistered().run(), map  // ignore unknown options
	);
	po::notify(map);

	if (map.count("showProbeDiscretization"))
	{
		cout << "openEMS - showing probe discretization information" << endl;
		m_showProbeDiscretization = true;
	}
	if (map.count("nativeFieldDumps"))
	{
		cout << "openEMS - dumping all fields using the native field components" << endl;
		m_nativeFieldDumps = true;
	}
	if (map.count("verbose"))
	{
		m_VerboseLevel = map["verbose"].as<unsigned int>();
	}

	if (m_VerboseLevel > 0) {
		cout << "openEMS - verbose level " << m_VerboseLevel << endl;
	}
}
