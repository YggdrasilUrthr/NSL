#pragma once 

#include <fstream>
#include <iostream>
#include <vector>

class csvwriter
{

private:

    std::string m_filename;
    std::ofstream m_fs;

public:

    csvwriter(const std::string filename);
    ~csvwriter();

    template<typename T> void write_csv_line(const std::vector<T> csv_line);
    void change_file(const std::string filename);

};

template<typename T> void csvwriter::write_csv_line(const std::vector<T> csv_line) {

    if(!m_fs.is_open()) {

        return;

    }

    for (T item : csv_line) {

        m_fs << item << ",";

    }

    m_fs.seekp(-1, std::ios_base::end);
    m_fs << std::endl;


}