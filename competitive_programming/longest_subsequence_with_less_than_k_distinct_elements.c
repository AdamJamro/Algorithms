//
// Created by adame on 3/24/2025.
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ALPHABET_SIZE 26


// -1 indicates error
// 0 indicates success
int map_add(int* map, const char key, int amount) {
    const int index = key - 'a';
    if (index >= ALPHABET_SIZE || key < 0) {
        return -1;
    }

    map[index] += amount;

    return 0;
}

// -1 indicates error
// x >= 0 is the value of the key value pair
int map_at(const int* map, const char key) {
    const int index = key - 'a';
    if (index >= ALPHABET_SIZE || key < 0) {
        return -1;
    }

    return map[index];
}

size_t max(const size_t a, const size_t b) {
    if (a > b) {
        return a;
    }
    return b;
}


// assumes sequence is lowercase letters a-z
size_t find_longest_subsequence(const char* sequence, const int k) {
    if (k <= 0) {
        return 0;
    }
    const size_t size = strlen(sequence);

    if (size == 0) {
        return 0;
    }
    if (size == 1) {
        return 1;
    }

    // simulate a hash map
    int* frequency_map = (int*) malloc(ALPHABET_SIZE * sizeof(int));
    for (int i = 0; i<ALPHABET_SIZE; i++) {
        frequency_map[i] = 0;
    }

    // here the following is true:
    // k is at least 1
    // sequence has at least two letters

    size_t left = 0;
    size_t right = 0;
    size_t max_length = 0; // imagine starting with range from 0 (inclusive) to 0 (exclusive) that has length of 0
    int variety_count = 1;
    map_add(frequency_map, sequence[0], 1);

    size_t index = right;
    while (index < size - 1) {
        while (variety_count <= k && right < size - 1) {
            right++;
            const size_t frequency = map_at(frequency_map, sequence[right]);
            if (frequency == 0) {
                variety_count++;
            }
            map_add(frequency_map, sequence[right], 1);
        }
        max_length = max(max_length, right - left); // right pointer is exclusive

        while (variety_count > k && left < right) {
            const size_t frequency = map_at(frequency_map, sequence[left]);
            map_add(frequency_map, sequence[left], -1);
            if (frequency == 1) {
                variety_count--;
            }
            left++;
        }

        // here we could either have:
        // 1) [left,right) point to a valid subsequence, and variety count equlas k => go to next iteration
        // 2) right == size - 1, and we have a valid solution including last element on [left, right + 1) => update max and break the loop
        // 3) left == right + 1 so [left,right) is an empty array that could only happen when k = 1 => just continue
        // 4) both [left,right) is an empty array and right == size - 1 =>
        // it's probably impossible after we handled sequences with size 1 separately, nevertheless just update the max and break the loop

        // 2) and 4) cases
        if (right == size - 1) {
            max_length = max(max_length, right - left + 1);
            break;
        }

        index = right;
    }

    // don't return inside the loop lest we could forget to free up the resources
    free(frequency_map);

    return max_length;
}

// success == 0; failure == -1;
int test(const char* sequence, const int letter_bound, const size_t expected_value) {
    if (expected_value != find_longest_subsequence(sequence, letter_bound)) {
        return -1;
    }
    return 0;
}

int main() {
    // char *foo = "aabbcd";
    // int k = 1;
    // printf("longest subsequence with at most %d different letters on a string %s is: %lu\n", k, foo, find_longest_subsequence(foo, k));

    int sum = 0;
    sum += test("ababbacabbaababbbaba", 2, 13);
    sum += test("aabbcd", 1, 2);
    sum += test("xyzxyzxyzxyz", 2, 2);
    sum += test("abbbaabcbaab", 3, 12);
    sum += test("abcdefghijklmno", 14, 14);
    sum += test("", 1, 0);
    sum += test("a", 1, 1);
    sum += test("ab", 0, 0);
    sum += test("ab", 1, 1);

    if (sum != 0) {
        printf("%d tests failed to pass\n", sum);
        return -1;
    }
    printf("All tests pass\n");
    return 0;
}

