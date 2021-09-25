//
// Created by 159-egi on 8/29/17.
//

#include "core/cell/side.h"

template <Side S>
inline constexpr const char* side_name_template() noexcept;

template <>
inline constexpr const char* side_name_template<Side::LEFT>() noexcept { return "left"; }

template <>
inline constexpr const char* side_name_template<Side::RIGHT>() noexcept { return "right"; }

template <>
inline constexpr const char* side_name_template<Side::TOP>() noexcept { return "top"; }

template <>
inline constexpr const char* side_name_template<Side::BOTTOM>() noexcept { return "bottom"; }

template <>
inline constexpr const char* side_name_template<Side::FRONT>() noexcept { return "front"; }

template <>
inline constexpr const char* side_name_template<Side::BACK>() noexcept { return "back"; }

const array<const char*, 6> side_to_string_vector = {
        side_name_template<to_side(0)>(),
        side_name_template<to_side(1)>(),
        side_name_template<to_side(2)>(),
        side_name_template<to_side(3)>(),
        side_name_template<to_side(4)>(),
        side_name_template<to_side(5)>()
};

template <Side S>
inline constexpr Side opposite_side_template() noexcept ;

template <>
inline constexpr Side opposite_side_template<Side::LEFT>() noexcept { return Side::RIGHT; }

template <>
inline constexpr Side opposite_side_template<Side::RIGHT>() noexcept { return Side::LEFT; }

template <>
inline constexpr Side opposite_side_template<Side::BOTTOM>() noexcept { return Side::TOP; }

template <>
inline constexpr Side opposite_side_template<Side::TOP>() noexcept { return Side::BOTTOM; }

template <>
inline constexpr Side opposite_side_template<Side::FRONT>() noexcept { return Side::BACK; }

template <>
inline constexpr Side opposite_side_template<Side::BACK>() noexcept { return Side::FRONT; }

const array<Side, 6> opposite_side_vector = {
        opposite_side_template<to_side(0)>(),
        opposite_side_template<to_side(1)>(),
        opposite_side_template<to_side(2)>(),
        opposite_side_template<to_side(3)>(),
        opposite_side_template<to_side(4)>(),
        opposite_side_template<to_side(5)>()
};

const static_vector<Side, 6> TWO_DIMENSIONAL_SIDES = {
        Side::LEFT,
        Side::BOTTOM,
        Side::RIGHT,
        Side::TOP
};


const static_vector<Side, 6> TWO_DIMENSIONAL_X_SIDES = {
        Side::LEFT,
        Side::RIGHT
};


const static_vector<Side, 6> TWO_DIMENSIONAL_Y_SIDES = {
        Side::BOTTOM,
        Side::TOP
};

const static_vector<Side, 6> THREE_DIMENSIONAL_SIDES = {
        Side::LEFT,
        Side::BOTTOM,
        Side::RIGHT,
        Side::TOP,
        Side::BACK,
        Side::FRONT
};