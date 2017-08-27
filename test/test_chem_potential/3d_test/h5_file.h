
#ifndef h5_file_h
#define h5_file_h

#include <string>
#include <vector>

#include "hdf5.h"

class H5File {
    private:
        hid_t m_file_id;

        template <typename T>
        hid_t getH5_Datatype();
        void create_group(std::string path);
        void parse (std::string attr_name, std::string &path, std::string &name);

    public:

        struct Error;

        H5File();

        void open(std::string filename, std::string mode);

        template <typename T>
        void read_dataset(std::string dataset_name, T * dataset);

        template <typename T>
        void write_dataset(std::string dataset_name, T * dataset, int * dims, int ndims);

        void get_ndims(std::string dataset_name, int &ndims);
        void get_dims(std::string dataset_name, int * dims);

        void list(std::string path, std::vector<std::string> &list);
        template <typename T>
        void set_attribute(std::string attr_name, T  attr_value);
        template <typename T>
        void get_attribute(std::string attr_name, T &attr_value);

        void close();
};

struct H5File :: Error : public std::exception {
    std::string msg;
    Error(std::string str) : msg(str) {}
    const char * what () const throw () {return msg.c_str();}
    ~Error() throw () {}
};

template <typename T>
hid_t H5File :: getH5_Datatype() {}
template <>
hid_t H5File :: getH5_Datatype<double>() { return H5T_NATIVE_DOUBLE; }
template <>
hid_t H5File :: getH5_Datatype<float>() { return H5T_NATIVE_FLOAT; }
template <>
hid_t H5File :: getH5_Datatype<int>() { return H5T_NATIVE_INT; }
template <>
hid_t H5File :: getH5_Datatype<unsigned int>() { return H5T_NATIVE_UINT; }

H5File :: H5File () {
    m_file_id = 0;
}

void H5File :: create_group(std::string path)
{
    int start = 1;
    size_t pos = 0;
    if (path[0]!='/') path = "/" + path;

    while ((pos = path.find("/", start)) != std::string::npos)
    {
        std::string substr = path.substr(0, pos);
        htri_t group_exists = H5Lexists(m_file_id, substr.c_str(), H5P_DEFAULT);
        if (!group_exists) {
            hid_t grp_id = H5Gcreate(m_file_id, substr.c_str() , H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Gclose(grp_id);
        }
        start = pos+1;
    }
}

void H5File :: parse (std::string attr_name, std::string &path, std::string &name)
{
    size_t start = 0;
    size_t end = 0;

    while (end != std::string::npos)
    {
        end = attr_name.find("/", start);
        path = attr_name.substr(0,start);
        name = attr_name.substr(start,end-start);
        start = end + 1;
    }
}

void H5File :: open(std::string filename, std::string mode)
{
    unsigned read, write, append;

    read = (mode == "r");
    write = (mode == "w");
    append = (mode == "a");

    // check that one is true
    if ( (read || write || append) != 1) 
        throw Error("Open file with \"r\", \"w\", or \"a\" mode");

    // check that a file is not already open
    if (m_file_id != 0) 
        throw Error("File is already open");

    if ( write ) {
        m_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else if ( read || append ) {
        FILE * file_exists = fopen(filename.c_str(), "r");

        if (file_exists == NULL) throw Error("File " + filename + " does not exists");
        else fclose(file_exists);

        if ( read )
            m_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if ( append )
            m_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

}

template <typename T>
void H5File :: write_dataset(std::string dataset_name, T * dataset, int * dims, int ndims)
{

    hid_t data_id;
    herr_t error;

    create_group(dataset_name);
    htri_t dataset_exists = H5Lexists(m_file_id, dataset_name.c_str(), H5P_DEFAULT);

    if ( dataset_exists == 0 ) {

        hsize_t h5_ndims = (hsize_t) ndims;
        hsize_t h5_dims[ndims];
        hid_t dcpl_id, space_id;

        // convert datatypes to be compatible with hdf5
        for (int i=0; i<ndims; i++) h5_dims[i] = (hid_t) dims[i];

        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl_id, h5_ndims, h5_dims);
        H5Pset_shuffle(dcpl_id);
        H5Pset_deflate(dcpl_id, 1);

        space_id = H5Screate_simple(h5_ndims, h5_dims, NULL);
        data_id = H5Dcreate(m_file_id,
                            dataset_name.c_str(),
                            getH5_Datatype<T>(),
                            space_id,
                            H5P_DEFAULT,
                            dcpl_id,
                            H5P_DEFAULT);

        H5Pclose(dcpl_id);
        H5Sclose(space_id);
        
    } else {
        data_id = H5Dopen(m_file_id, dataset_name.c_str(), H5P_DEFAULT);
    }

    error = H5Dwrite(data_id,
                     getH5_Datatype<T>(),
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     dataset);

    H5Dclose(data_id);

    if ( dataset_exists )
        throw Error("Trying to overwrite dataset " + dataset_name);
    if ( error < 0 )
        throw Error("Error writing dataset " + dataset_name);

}

void H5File :: get_ndims(std::string dataset_name, int &ndims)
{
    hid_t data_id, space_id;

    data_id = H5Dopen(m_file_id, dataset_name.c_str(), H5P_DEFAULT);
    space_id = H5Dget_space(data_id);
    ndims = H5Sget_simple_extent_ndims(space_id);

    H5Sclose(space_id);
    H5Dclose(data_id);

}

void H5File :: get_dims(std::string dataset_name, int * dims)
{
    int ndims;
    int max_dims=10;
    hid_t data_id, space_id;
    hsize_t h5_dims[max_dims];

    data_id = H5Dopen(m_file_id, dataset_name.c_str(), H5P_DEFAULT);
    space_id = H5Dget_space(data_id);
    ndims = H5Sget_simple_extent_ndims(space_id);
    if (ndims > max_dims) 
        throw Error("Error with number of dimensions ");

    H5Sget_simple_extent_dims(space_id, h5_dims, NULL);
    for (int i=0; i<ndims; i++) dims[i] = (int) h5_dims[i];

    H5Sclose(space_id);
    H5Dclose(data_id);

}

template <typename T>
void H5File :: read_dataset(std::string dataset_name, T * dataset)
{
    hid_t data_id;
    herr_t error;

    data_id = H5Dopen(m_file_id, dataset_name.c_str(), H5P_DEFAULT);
    error = H5Dread(data_id, getH5_Datatype<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset);
    if (error < 0) 
        throw Error("Error reading dataset " + dataset_name);

    H5Dclose(data_id);
}

void H5File :: list(std::string path, std::vector<std::string> &list)
{

    /**
    @param path list all the dataset in this path
    @param list modified to contain the names of the datasets
    **/

    hid_t group_id;
    H5G_info_t group_info;
    int nlinks;

    // need to implement group existence check

    group_id = H5Gopen(m_file_id, path.c_str(), H5P_DEFAULT);
    H5Gget_info(group_id, &group_info);
    nlinks = group_info.nlinks;

    list.resize(nlinks);

    for (int i=0; i<nlinks; i++)
    {
        char name[256];
        H5Gget_objname_by_idx(group_id, i, name, 256);
        list[i] = name;
    }


    H5Gclose(group_id);

}

template <typename T>
void H5File :: set_attribute(std::string attr_name, T  attr_value)
{
    /**
    * @param attr_name Name/address of the attribute
    * @param attr_value Value of the attribute
    **/

    /**
    The attr_name must have atleast one backslash.
    For example, setting an attribute named my_attr in the root group
    should be specified as "/my_attr"
    */

    /** \error ERROR 1: attr_name not sufficient, ensure that it begins with a backslash */
    /** \error ERROR 2: attr_name not sufficient, ensure that it does not end in a backslash */

    std::string path;
    std::string name;
    parse(attr_name, path, name);

    if (path.size() == 0 || name.size() == 0)
        throw Error("Error setting attribute " + attr_name);

    create_group(attr_name);
    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate_by_name(m_file_id, 
                              path.c_str(),
                              name.c_str(), 
                              getH5_Datatype<T>(), 
                              space_id, 
                              H5P_DEFAULT, 
                              H5P_DEFAULT,
                              H5P_DEFAULT);

    H5Awrite(attr_id, getH5_Datatype<T>(), &attr_value);
    H5Aclose(attr_id);
    H5Sclose(space_id);

}

template <typename T>
void H5File :: get_attribute(std::string attr_name, T &attr_value)
{
    /**
    * @param attr_name Name/address of the attribute
    * @param attr_value Value of the attribute
    **/

    /** \error ERROR 1: attribute does not exist */

    hid_t attr_id;
    std::string path;
    std::string name;

    parse(attr_name, path, name);

    if (H5Aexists_by_name(m_file_id, path.c_str(), name.c_str(), H5P_DEFAULT) <= 0) 
        throw Error("Attempting to get attribute which does not exist: " + attr_name);

    attr_id = H5Aopen_by_name(m_file_id, path.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    H5Aread(attr_id, getH5_Datatype<T>(), &attr_value);

    H5Aclose(attr_id);

}

void H5File :: close()
{
    if ( m_file_id == 0 ) {
        throw Error("Attempting to close file that is not open");
    } else {
        H5Fclose(m_file_id);
        m_file_id = 0;
    }
}

#endif
