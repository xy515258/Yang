
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>

class ParameterFile {

    private:

        std::map<std::string, std::string> params;

        void removeWhitespace(std::string & str) const;
        bool checkVectorSyntax(const std::string & vec_string) const;
        std::vector<std::string> parse(const std::string & str, const char delim) const;

        template<typename T>
        T convertValueType(const std::string str) const;

    public:

        struct Error;

        void readParameters(const std::string & filename);

        std::map<std::string, std::string> getParameters() const;

        template<typename T>
        void unpack( const std::string & name, T & parameter  ) const;

        template<typename T>
        void unpack( const std::string & name, std::vector<T> & vec  ) const;

};

struct ParameterFile :: Error : public std::exception {
    std::string msg;
    Error(std::string str) : msg(str) {}
    const char * what () const throw () {return msg.c_str();}
    ~Error() throw () {}
};

template<typename T>
T ParameterFile :: convertValueType(const std::string str) const
{
    T value;
    std::stringstream ss(str);
    ss >> value;
    if (ss.fail()) throw Error("Invalid input: " + str);
    return value;
}

template<typename T>
void ParameterFile :: unpack( const std::string & name, T & parameter  ) const
{
    std::map<std::string,std::string>::const_iterator it = params.find(name);
    if (it == params.end()) throw Error("Parameter not found: " + name);

    parameter = convertValueType<T>(it->second);
}

template<typename T>
void ParameterFile :: unpack( const std::string & name, std::vector<T> & vec  ) const
{
    // check that the parameter is found
    std::map<std::string,std::string>::const_iterator it = params.find(name);
    if (it == params.end()) throw Error("Parameter not found: " + name);

    // check the syntax
    std::string vec_string = it->second;
    bool isVector = checkVectorSyntax(vec_string);
    if (isVector == false) throw Error("Improper vector syntax: " + vec_string);

    // parse vector elements
    std::string no_braces = vec_string.substr(1,vec_string.length()-2);
    std::vector<std::string> list = parse(no_braces, ',');

    for (unsigned i=0; i<list.size(); i++)
    {
        T value = convertValueType<T>(list[i]);
        vec.push_back(value);
    }
}

#endif
