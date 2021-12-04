#pragma once

#include <string>
#include <boost/algorithm/string.hpp>


/// @brief Флаги граничных условий, устанавливаются с двух сторон от граней.
enum class FaceFlag : int {
    UNDEFINED = 0,   ///< Не определено.
    ORDER = 1,       ///< Обычное соседство (внутри сетки).
    WALL = 2,        ///< Непроницаемая стенка.
    ZOE = 3,         ///< Отображение данных запрашиваемой ячейки (простой снос).
    PERIODIC = 4,    ///< Периодичность.
    FUNCTION = 5,    ///< Функция задается пользователем.
    INFLOW = 6,      ///< Втекание жидкости.
    OUTFLOW = 7,     ///< Вытекание жидкости.
};

/// @brief Преобразует флаг в целочисленный тип.
inline int to_int(FaceFlag flag) {
    return static_cast<int>(flag);
}

/// @brief Преобразует целое число во флаг.
inline FaceFlag to_face_flag(int flag) {
    return static_cast<FaceFlag>(flag);
}

/// @brief Преобразует флаг в строку.
inline std::string to_string(FaceFlag flag) {
    switch (flag) {
        case FaceFlag::ORDER :
            return "ORDER";
        case FaceFlag::WALL:
            return "WALL";
        case FaceFlag::ZOE:
            return "ZOE";
        case FaceFlag::FUNCTION:
            return "FUNCTION";
        case FaceFlag::PERIODIC:
            return "PERIODIC";
        case FaceFlag::UNDEFINED:
            return "UNDEFINED";
        case FaceFlag::INFLOW:
            return "INFLOW";
        case FaceFlag::OUTFLOW:
            return "OUTFLOW";
        default:
            throw std::runtime_error("Unknown FaceFlag");
    }
}

/// @brief Преобразует строку во флаг.
/// В таком виде необходимо задавать граничные условия в
/// конфигурационном файле.
inline FaceFlag to_face_flag(std::string str) {
    boost::algorithm::to_upper(str);

    if (str == "ORDER") {
        return FaceFlag::ORDER;
    } else if (str == "WALL") {
        return FaceFlag::WALL;
    } else if (str == "ZOE") {
        return FaceFlag::ZOE;
    } else if (str == "FUNCTION") {
        return FaceFlag::FUNCTION;
    } else if (str == "PERIODIC") {
        return FaceFlag::PERIODIC;
    } else if (str == "INFLOW") {
        return FaceFlag::INFLOW;
    } else if (str == "OUTFLOW") {
        return FaceFlag::OUTFLOW;
    } else {
        throw std::runtime_error("Error: Can't read FaceFlag from string " + str);
    }
}
