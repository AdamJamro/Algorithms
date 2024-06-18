#include <iostream>
#include <vector>

using std::string;
using std::vector;

vector<unsigned int> KMP(const string&, const string&);
unsigned int* compute_kmp_prefix_function(const string&);

// returns all starting positions where pattern occurs
vector<unsigned int> KMP(const string& text, const string& pattern) {
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
    offsets.reserve(text_len);

    auto* prefix = compute_kmp_prefix_function(pattern);
    unsigned int letters_matched = 0;
    for (int i = 0; i < text_len; i++) {
        while (letters_matched > 0 and pattern[letters_matched] != text[i]) {
            letters_matched = prefix[letters_matched - 1];
        }

        if (pattern[letters_matched] == text[i]){
            letters_matched++;
        }

        if (letters_matched == pattern_len) {
            offsets.push_back(i - pattern_len + 1);
            letters_matched = prefix[letters_matched - 1];
        }
    }

    return offsets;
}

unsigned int* compute_kmp_prefix_function(const string& pattern) {
    unsigned int len = pattern.length();
    if (len < 1) throw std::invalid_argument("empty string");

    auto* prefix = new unsigned int[len];
    prefix[0] = 0;
    unsigned int letters_matched = 0;
    for (int i = 1; i < len; i++) {
        while (letters_matched > 0 and pattern[letters_matched] != pattern[i]) {
            letters_matched = prefix[letters_matched - 1];
        }

        if (pattern[letters_matched] == pattern[i]){
            letters_matched++;
        }

        prefix[i] = letters_matched;
    }
    return prefix;
}

void print_string(int from, int to, const string& s) {
    for (int i = from; i < to; i++) {
        std::cout << s[i];
    }
}

int main() {
    string tex = "00110011000101001110001110101";
    string pattern = "0100";
    for (const auto& index : KMP(tex, pattern)) {
        if (index < 0) {
            std::cout << "no matches";
        }
        else {
            std::cout<<index << ": ";
            print_string(index, index + pattern.length(), tex);
            std::cout<<"\n";
        }
    }

    return 0;

}
