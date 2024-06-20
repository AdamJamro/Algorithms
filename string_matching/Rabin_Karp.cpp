#include <iostream>
#include <vector>
#include <fstream>

using std::string;
using std::vector;

vector<unsigned int> rabin_karp(const string&, const string&, unsigned int, unsigned int);
long int pow(long int, long int);
inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}

long int pow(long int base, long int exponent) {
    int x = base;
    int res = 1;
    while (exponent > 0){
        if (exponent & 1) {
            res = res * x;
        }
        x = x * x;
        exponent = exponent >> 1;
    }
    return res;
}

vector<unsigned int> rabin_karp(const string& text, const string& pattern, unsigned int letter_size, const unsigned int modular = 2011) {
    unsigned int text_len = text.length();
    unsigned int pattern_len = pattern.length();

    vector<unsigned int> offsets{};
    if (pattern_len == 0){
        offsets.push_back(0);
        return offsets;
    }
    if (text_len < pattern_len) {
        offsets.push_back(-1);
        return offsets;
    }
    offsets.reserve(int(text_len / pattern_len));

    unsigned int h = pow(letter_size, pattern_len - 1) % modular;
    unsigned int integer_pattern = 0; // integer_pattern := pattern converted into integer
    unsigned int t_0 = 0; // t_s := substring of text[s, s+1, ..., s + pattern_len - 1] converted to integer

    // preprocessing (Horner's rule)
    for (int i = 0; i < pattern_len; i++) {
        integer_pattern = (letter_size * integer_pattern + pattern[i]) % modular;
        t_0 = (letter_size * t_0 + text[i]) % modular;
    }
    // matching
    unsigned int t_offset = t_0;
    for (auto offset = 0; offset < text_len - pattern_len + 1; offset++) {
        if (integer_pattern == t_offset) {
            if (pattern == text.substr(offset,  pattern_len)) {
                offsets.push_back(offset);
            }
        }
        if (offset < text_len - pattern_len + 1) {
            t_offset = positive_modulo(letter_size * (t_offset - text[offset] * h) + text[offset + pattern_len], modular);
        }
    }

    if (offsets.empty()) {
        offsets.push_back(-1);
    }
    return offsets;
}

//void print_string(int from, int to, const string& s) {
//    for (int i = from; i < to; i++) {
//        std::cout << s[i];
//    }
//}

int main(int argc, char *argv[]) {
    string error_line = "Run this file with arguments KMP [pattern] [file_name]\n";
    if (argc != 3) {
        std::cout << "this program accepts exactly 2 arguments!\n" << error_line;
        return -1;
    }
    string pattern = argv[1];
    string file_name = argv[2];
    if (pattern.empty() or file_name.empty()) {
        std::cout << "invalid arguments!\n" << error_line;
        return -1;
    }

    string line;
    long int line_count = 0;
    std::ifstream file_stream (file_name);

    if (file_stream.is_open()) {

        while ( getline (file_stream, line) )
        {
            line_count++;
            vector<unsigned int> found_matches = rabin_karp(line, pattern, pow(2, sizeof(char) * 8));

            if (found_matches[0] != -1) {
                std::cout << "Line " << line_count << " at: ";

                for (const auto& offset : found_matches) {
                    std::cout << offset << " ";
                }
                std::cout << "\n";
            }
        }
        file_stream.close();
    }

    else std::cout << "Unable to open file";

    return 0;
}
