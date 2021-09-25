//
// Created by 159-mrv on 4/16/18.
//

#ifndef INFINIMESH_DATA_HOLDER_H
#define INFINIMESH_DATA_HOLDER_H

#include <allstd.h>

/**
 * @brief Абстрактный контейнер для хранения данных типа double
 */
class DataHolder {
public:
    using Ptr = shared_ptr<DataHolder>;
    using Ref = const shared_ptr<DataHolder> &;

    explicit DataHolder(size_t size);

    static DataHolder::Ptr create(size_t size) {
        return make_shared<DataHolder>(size);
    }

    /** @brief Размер данных в double */
    size_t size() const;

    /** @brief Расширить буффер данных */
    void resize(size_t new_size);

    /**
     * @brief Перенести данные в буффер
     * @param output буффер для сохранения данных
     */
    void serialize(double *output) const;

    /**
     * @brief Записать данные из буффера
     * @param input буффер-источник данных
     */
    void deserialize(const double *input);

    /** @return Указатель на буффер данных */
    double* buffer();

    /** @return Хэш данных */
    size_t hash() const;

protected:
    vector<double> m_data;
};

#endif //INFINIMESH_DATA_HOLDER_H
