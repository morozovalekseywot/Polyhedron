//
// Created by 159-mrv on 10/27/17.
//

#ifndef AMR_HDF_FILE_H
#define AMR_HDF_FILE_H

#include <allstd.h>

#include <H5Ipublic.h>

namespace HDF {

    class Object {
    public:
        virtual void close() = 0;
    protected:
        hid_t m_id;
        explicit Object(hid_t id = -1) : m_id(id) {}
    };

    class Group : public Object {
    public:
        Group(Group* parent, const string& name);

        Group add_group(const string& name);
        void add_attribute(const string& name, hid_t type, void* data);
        void add_table(const string& name, hid_t type, const vector<hsize_t>& sizes);

        Group subgroup(const string& name);

        void write(const string& table_name, hid_t type, void * data);

        void close() override;

    protected:
        explicit Group(hid_t id = -1) : Object(id) {};

    };


    class File : public Group {
    public:
        File() = default;
        explicit File(const string& filename);
        void close() override;

    public:
        static void create(const string& filename);
    };

}
#endif //AMR_HDF_FILE_H
