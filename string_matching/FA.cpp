#include <iostream>
#include <vector>
#include <fstream>

using std::string;
using std::vector;

vector<unsigned int> finite_automaton_matcher(const string&, const string&, size_t);
vector<unsigned int> finite_automaton_matcher(const string&, const string&, size_t, const unsigned int*);
unsigned int* compute_transition_function(const string&, size_t);

inline size_t make_2d(size_t a, size_t b, size_t alphabet_size=256) {
    return a * alphabet_size + b;
}

// function allocates memory and never deletes it
unsigned int* compute_transition_function(const string& pattern, const size_t alphabet_size=256) {
    const size_t pattern_len = pattern.length();

    auto* prefix = new unsigned int[pattern_len];
    prefix[0] = 0;
    unsigned int letters_matched = 0;
    for (int i = 1; i < pattern_len; i++) {
        while (letters_matched > 0 and pattern[letters_matched] != pattern[i]) {
            letters_matched = prefix[letters_matched - 1];
        }

        if (pattern[letters_matched] == pattern[i]){
            letters_matched++;
        }

        prefix[i] = letters_matched;
    }

    auto* transition_fun = new unsigned int[(pattern_len + 1) * alphabet_size];

    for (int letter = 0; letter < alphabet_size; letter++) {
        if (letter == pattern[0]) {
            transition_fun[make_2d(0, letter, alphabet_size)] = 1;
        }
        else {
            transition_fun[make_2d(0, letter, alphabet_size)] = 0;
        }
        for (int i = 1; i <= pattern_len; i++) {
            const int& state = i;
            letters_matched = state; // imagine we're already in the ith state with i letters matched

            // try to see if current prefix matches against current letter
            while (letters_matched > 0 and pattern[letters_matched] != letter) {
                letters_matched = prefix[letters_matched - 1]; // precomputed prefix gives returns only next valid shifts
            }

            if (pattern[letters_matched] == letter) {
                letters_matched++;
            }

            // from current state with current letter we go into ith state
            // (we found input postfix of length i that matches ith prefix of pattern)
            auto& next_state = transition_fun[make_2d(state, letter, alphabet_size)];
            if (pattern_len == state and letters_matched == pattern_len) {
                next_state = 0;
            }
            else {
                next_state = letters_matched;
            }
        }
    }

    delete[] prefix;
    return transition_fun;
}

vector<unsigned int> finite_automaton_matcher(const string& text, const string& pattern, const size_t alphabet_size=256) {
    const size_t text_len = text.length();
    const size_t pattern_len = pattern.length();

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


    // creating the automaton
    auto* transition_fun = compute_transition_function(pattern, alphabet_size);


    // simulation
    unsigned int q = 0;
    for (int i = 0; i < text_len; i++) {
        q = transition_fun[make_2d(q, text[i], alphabet_size)];

        if (q == pattern_len) {
            offsets.push_back(i - pattern_len + 1);
        }
    }

    if (offsets.empty()) {
        offsets.push_back(-1);
    }

    delete transition_fun;
    return offsets;
}

vector<unsigned int> finite_automaton_matcher(const string& text, const string& pattern, const size_t alphabet_size, const unsigned int* transition_fun) {
    const size_t text_len = text.length();
    const size_t pattern_len = pattern.length();

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

    // simulation
    unsigned int q = 0;
    for (int i = 0; i < text_len; i++) {
        q = transition_fun[make_2d(q, text[i], alphabet_size)];

        if (q == pattern_len) {
            offsets.push_back(i - pattern_len + 1);
        }
    }

    if (offsets.empty()) {
        offsets.push_back(-1);
    }

    return offsets;
}

int main(int argc, char *argv[]) {
    string error_line = "Run this file with arguments FA [pattern] [file_name]\n";
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

    size_t alphabet_size = 256;
    // precompute the automaton
    auto* transition_function = compute_transition_function(pattern, 256);


    string line;
    long unsigned int line_count = 0;
    std::ifstream file_stream (file_name);
    if (file_stream.is_open()) {

        while ( getline (file_stream, line) )
        {
            line_count++;
            vector<unsigned int> found_matches = finite_automaton_matcher(line, pattern, alphabet_size, transition_function);

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
    delete[] transition_function;
    return 0;
}
