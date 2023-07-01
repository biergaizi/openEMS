/*
 * Copyright (C) 2023 Yifeng Li <tomli@tomli.me>
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

#include <vector>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>

#ifndef TILING_H
#define TILING_H

typedef std::pair<int, int>	Range;
typedef std::vector<Range>	Block;

enum {
	TILES_RECTANGULAR,
	TILES_PARALLELOGRAM,
	TILES_DIAMOND,
};

struct Tiles {
	/* TILES_PARALLELOGRAM or TILES_DIAMOND */
	int type;

	/*
	 * TILES_PARALLELOGRAM has one phases,
	 * TILES_DIAMOND has two.
	 */
	int phases;

	/* 
	 * First index: phase.
	 * TILES_PARALLELOGRAM has one phase.
	 * TILES_DIAMOND has mountain and valley phases.
	 *
	 * Second index: block.
	 * Depending on the tiling parameters, the width
	 * of a block varies.
	 *
	 * Third index: timestep.
	 * Depending on the tiling parameters, the number
	 * of timesteps within one block varies.
	 */
	std::vector<std::vector<Block>> array;
};

Tiles computeParallelogramTiles1D(
	int totalWidth,
	int blkWidth,
	int blkTimesteps
);
Tiles computeDiamondTiles1D(
	int totalWidth,
	int blkWidth,
	int blkTimesteps
);

struct Range3D
{
	int timestep;
	int voltageStart[3];
	int voltageStop[3];
	int currentStart[3];
	int currentStop[3];
};
typedef std::vector<Range3D>	Tiles3D;

Tiles computeRectangularTilesNoDeps1D(int totalWidth, int blkWidth, int blkHalfTimesteps);
Tiles computeRectangularTiles1D(int totalWidth, int blkWidth, int blkHalfTimesteps);
std::vector<std::vector<Tiles3D>> computeRectangularTiles3D(
	int totalWidth[3],
	int blkWidth[3],
	int numThreads
);
Tiles3D combineTilesTo3D(Tiles tilesX, Tiles tilesY, Tiles tilesZ, int blkHalfTimesteps);
void visualizeTiles(Tiles tiles, int totalWidth, int blkTimesteps);
std::vector<std::vector<Tiles3D>> combineTilesTo3D(
	Tiles tilesX, Tiles tilesY, Tiles tilesZ,
	int blkHalfTimesteps,
	int numThreads
);

#endif
