//
// Created by 159-mrv on 7/12/18.
//

#include <core/cell/data_holder.h>


DataHolder::DataHolder(size_t size)
    : m_data(vector<double>(size, 0.0)) {
}

size_t DataHolder::size() const {
    return m_data.size();
}

void DataHolder::resize(size_t new_size) {
    m_data.resize(new_size);
}

void DataHolder::serialize(double *output) const {
    memcpy(output, m_data.data(), m_data.size() * sizeof(double));
}

void DataHolder::deserialize(const double *input) {
    memcpy(m_data.data(), input, m_data.size() * sizeof(double));
}

double * DataHolder::buffer() {
    return m_data.data();
}

size_t DataHolder::hash() const {
    size_t result = 0;
    for(double d: m_data) {
         result = result * 31 + std::hash<double>()(d);
    }
    return result;
}