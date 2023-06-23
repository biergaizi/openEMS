#include <iostream>
#include <tuple>
#include <unordered_map>
#include <boost/container_hash/extensions.hpp>

struct Block {
        int start;
        int end;
};

// std::tuple<order, start[3], stop[3]>
using TileKey = std::tuple<int, int*, int*>;

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

//using Coords3 = std::array<unsigned int, 3>;
using TileMap = std::unordered_map<TileKey, std::vector<int>, TileKeyHash>;
