//
// Created by egi on 8/27/17.
//

#ifndef INFINIMESH_SIDE_H
#define INFINIMESH_SIDE_H

#include <allstd.h>

/// @brief Стороны ячейки
enum class Side : int
{
    LEFT = 0,
    BOTTOM = 1,
    RIGHT = 2,
    TOP = 3,
    FRONT = 4,
    BACK = 5,
    CENTER = 6
};

inline constexpr int to_int(Side side) noexcept {
    return static_cast<int>(side);
}

inline constexpr Side to_side(int idx) noexcept {
    return static_cast<Side>(idx);
}

extern const array<const char*, 6> side_to_string_vector;

inline const char* to_string(Side side) noexcept {
    return side_to_string_vector[to_int(side)];
}

extern const array<Side, 6> opposite_side_vector;

inline Side opposite_side(Side side) noexcept {
    return opposite_side_vector[to_int(side)];
}

extern const static_vector<Side, 6> TWO_DIMENSIONAL_SIDES;

extern const static_vector<Side, 6> THREE_DIMENSIONAL_SIDES;

extern const static_vector<Side, 6> TWO_DIMENSIONAL_X_SIDES;

extern const static_vector<Side, 6> TWO_DIMENSIONAL_Y_SIDES;

#endif // INFINIMESH_SIDE_H
