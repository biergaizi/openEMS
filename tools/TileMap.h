#include <iostream>
#include <tuple>
#include <unordered_map>
#include <boost/container_hash/extensions.hpp>

struct Block {
        int start;
        int end;
};

using TileKey = std::tuple<int, std::array<int, 3>, std::array<int, 3>>;

struct TileKeyHash {
	std::size_t operator()(const TileKey& key) const
	{
		std::tuple<int, int, int, int, int, int, int> t = {
			std::get<0>(key),
			std::get<1>(key)[0], std::get<1>(key)[1], std::get<1>(key)[2],
			std::get<2>(key)[0], std::get<2>(key)[1], std::get<2>(key)[2]
		};

		return boost::hash_value(t);
	}
};

inline TileKey GetTileKey(int order, int start[3], int end[3])
{
	std::array<int, 3> startArray = {
		start[0], start[1], start[2]
	};
	std::array<int, 3> stopArray = {
		end[0], end[1], end[2]
	};
	return std::make_tuple(order, startArray, stopArray);
}

//using Coords3 = std::array<unsigned int, 3>;
using TileMap = std::unordered_map<TileKey, std::vector<int>, TileKeyHash>;
