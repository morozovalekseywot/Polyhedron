//
// Created by 159-mrv on 10/27/17.
//

#include <iostream>
#include <hdf5.h>

#include "io/hdf/hdf_utils.h"

using namespace HDF;


void File::create(const string& filename) {
    hid_t hdf5 = H5Fcreate((filename + ".hdf5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(hdf5);
}

File::File(const string& filename) {
    m_id = H5Fopen((filename + ".hdf5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
}

Group::Group(Group *parent, const string& name) {
    m_id = H5Gcreate(parent->m_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

Group Group::add_group(const string& name) {
    return Group(this, name);
}

void Group::add_attribute(const string& name, hid_t type, void *data) {
    auto aid2 = H5Screate(H5S_SCALAR);
    auto attribute = H5Acreate2 (m_id, name.c_str(), type,
                                 aid2, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attribute, type, data);
    H5Sclose(aid2);
    H5Aclose(attribute);
}

void Group::add_table(const string& name, hid_t type, const vector<hsize_t>& sizes) {
    hid_t dataset, dataspace;
    dataspace = H5Screate_simple((int)sizes.size(), sizes.data(), NULL);
    dataset = H5Dcreate(m_id, name.c_str(), type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

Group Group::subgroup(const string& name) {
    return Group(H5Gopen(m_id, name.c_str(), H5P_DEFAULT));
}

void Group::write(const string& table_name, hid_t type, void * data) {
    auto data_set = H5Dopen(m_id, table_name.c_str(), H5P_DEFAULT);
    //auto data_space = H5Dget_space(data_set);
    //hsize_t dims[5];
    //H5Sget_simple_extent_dims(data_space, dims, NULL);
    H5Dwrite(data_set, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(data_set);
}

void Group::close() {
    H5Gclose(m_id);
}

void File::close() {
    H5Fclose(m_id);
}


