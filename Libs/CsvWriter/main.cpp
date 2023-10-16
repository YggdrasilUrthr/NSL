#include "csvwriter.h"

int main() {

    csvwriter Writer("test.csv");

    std::vector<uint32_t> test_line = {1, 2, 3, 4};
    std::vector<uint32_t> test_line_b = {5, 6, 7, 8};
    Writer.write_csv_line<uint32_t>(test_line);
    Writer.write_csv_line<uint32_t>(test_line_b);

    return 0;

}