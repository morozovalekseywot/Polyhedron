//
// Created by 159-mrv on 2/15/19.
//

#ifndef INFINIMESH_ORTH_DECOMPOSITION_H
#define INFINIMESH_ORTH_DECOMPOSITION_H

#include <core/generator/decomposition.h>

enum class GeometryType : int;

class SectorGenerator;
class RealSectorGenerator;
class RectangleGenerator;
class Sector3dGenerator;
class CuboidGenerator;

class OrthDecomposition : public Decomposition {

    enum class Direction : int { Ox = 0, Oy = 1, Oz = 2 };

public:

    explicit OrthDecomposition(const Configuration &config,
                               GeometryGenerator* geometry);

    int rank(const Vector3d &vertex) const final;


    void set_workload(double workload) final;

    void sync_data() final;

    void update() final;

    ~OrthDecomposition() final = default;


private:

    /// @brief Инициализирует поля m_dirs, m_px, m_py, m_pz
    void read_config(const Configuration &config);

    void init_barriers(GeometryGenerator* geometry);

    /// @brief Заполняет поля m_barriers для секторальных областей
    void sector_barriers_init(SectorGenerator* sector);

    /// @brief Заполняет поля m_barriers для секторальных областей
    void real_sector_barriers_init(RealSectorGenerator* sector);

    /// @brief Заполняет поля m_barriers для прямоугольной области
    void rectangle_barriers_init(RectangleGenerator* rectangle);

#if DIM3
    /// @brief Заполняет поля m_barriers для трехмерного сектора
    void sector3d_barriers_init(Sector3dGenerator* sector);

    /// @brief Заполняет поля m_barriers для кубоида
    void cuboid_barriers_init(CuboidGenerator* cuboid);
#endif

    void common_sector_barriers_init(double r0, double r1, double r2, double alpha);

    /// @brief Заполняет поля m_barriers
    /// @param get_coord - Отображение единичного квадрата/куба на
    /// область паралеллепипеда m_barriers
    void set_barriers(const std::function<double(Direction, double)> &get_coord);

    /// @brief Небольшое смещение барьеров для нарушения возможной симметрии
    void skew_barriers();

    /// @brief Возвращает координаты процесса по ранку
    array<int, 3> proc_coord(int rank) const;


    /// @brief Тип геометрии, необходим для начального расположения границ
    /// блоков и преобразования координат при смещении
    GeometryType m_type;

    /// @brief Порядок декомпозиции
    array<Direction, 3> m_dir;

    /// @brief Число процессов вдоль направления
    int m_px;
    vector<int> m_py;
    vector<vector<int>> m_pz;

    /// @brief Координаты текущего процесса
    int m_x, m_y, m_z;

    /// @brief Координаты перегородок
    vector<double> m_barriers_x;
    vector<vector<double>> m_barriers_y;
    vector<vector<vector<double>>> m_barriers_z;

    /// @brief Смещения перегородок
    vector<double> m_dx;
    vector<vector<double>> m_dy;
    vector<vector<vector<double>>> m_dz;

    /// @brief Нагрузка на процессы
    vector<vector<vector<double>>> m_workload;

    /// @brief Переход к координатам декомпозиции
    std::function<Vector3d(const Vector3d&)> m_to_local;

};

#endif //INFINIMESH_ORTH_DECOMPOSITION_H
