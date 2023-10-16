#include "csvwriter.h"

csvwriter::csvwriter(const std::string filename) {

    m_filename = filename;
    m_fs = std::ofstream(filename);

    if(!m_fs.is_open()) {

        std::cerr << "Unable to open file." << std::endl;        

    }

}

csvwriter::~csvwriter() {

    if(m_fs.is_open()) {

        m_fs.close();

    }
    

}

void csvwriter::change_file(const std::string filename) {

    if(m_fs.is_open()) {

        m_fs.close();

    }

    m_filename = filename;
    m_fs = std::ofstream(filename);

    if(!m_fs.is_open()) {

        std::cerr << "Unable to open file." << std::endl;        

    }

}