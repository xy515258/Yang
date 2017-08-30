
#include "parameter_file.h"

std::map<std::string, std::string> ParameterFile :: getParameters() const
{
    return params;
}


void ParameterFile :: removeWhitespace(std::string & str) const
{
    std::string::iterator end = str.end();
    end = std::remove(str.begin(),end,' ');
    end = std::remove(str.begin(),end,'\t');
    end = std::remove(str.begin(),end,'\n');
    end = std::remove(str.begin(),end,'\r');
    end = std::remove(str.begin(),end,'\v');
    str.erase(end, str.end());
}


bool ParameterFile :: checkVectorSyntax(const std::string & vec_string) const
{
    int n = vec_string.length();
    if (vec_string[0] != '{') return false;
    if (vec_string[n-1] != '}') return false;
    if (*(vec_string.begin()+1)==',') return false;
    if (*(vec_string.end()-2)==',') return false;
    for (int i=1; i<n-1; i++)
        if (vec_string[i]=='{' || vec_string[i]=='}') 
            return false;
    for (int i=1; i<n; i++)
        if (vec_string[i]==',' && vec_string[i-1]==',') 
            return false;
    return true;
}

std::vector<std::string> ParameterFile :: parse(const std::string & str, const char delim) const
{
    std::stringstream ss(str);
    std::vector<std::string> result;
    if (str.empty()) return result;
    while (ss.good())
    {
        std::string substr;
        std::getline( ss, substr, delim );
        result.push_back( substr );
    }
    return result;
}

void ParameterFile :: readParameters(const std::string & filename)
{
    /** 
    This function parses the input file.
    It reads a key-value pair from any line involving an "=" sign.
    All key-value pairs are stored as strings in params for use later.
    Comment lines begin with "#"
    */

    std::ifstream input(filename.c_str());
    std::string line;

    if (!input.good()) throw Error("File can't be found/opened: " + filename);

    while (std::getline(input, line))
    {
        size_t pos;
        std::string key, value;

        removeWhitespace(line);

        // skip empty lines
        if (line.empty()) continue;

        // skip comments
        if (line[0] == '#') continue;

        // parse line using the equal sign
        pos = line.find("=");
        key = line.substr(0, pos);
        value = line.substr(pos+1, line.length());

        params[key] = value;
    }

    input.close();
}
