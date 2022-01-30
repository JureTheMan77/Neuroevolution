//
// Created by jure on 4. 02. 21.
//

#ifndef NEUROEVOLUTION_UTIL_H
#define NEUROEVOLUTION_UTIL_H

#include <string>
#include <vector>
#include <iostream>
#include <random>

namespace util {

    inline std::random_device seeder;
    inline std::uniform_real_distribution<double> probabilityDistributionDouble(0, 1);
    inline std::uniform_int_distribution<unsigned int> probabilityDistributionInt(0, 1);

    /**
     * Split a string into tokens.
     * Source @see <a href="https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c">stackoverflow</a>
     * @param text Text to split
     * @param delimiter Delimiter to find in text
     * @return Vector of tokens
     */
    inline std::vector<std::string> split(std::string &text, std::string_view delimiter) {
        auto tokens = std::vector<std::string>();
        size_t pos;
        while ((pos = text.find(delimiter)) != std::string::npos) {
            tokens.push_back(text.substr(0, pos));
            text.erase(0, pos + delimiter.length());
        }
        tokens.push_back(text);
        return tokens;
    }

    /**
     * Split a string into tokens of type double.
     * Source @see <a href="https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c">stackoverflow</a>
     * @param text Text to split
     * @param delimiter Delimiter to find in text
     * @return Vector of double tokens
     */
    inline std::vector<double> splitDouble(std::string &text, std::string_view delimiter) {
        auto tokens = std::vector<double>();
        size_t pos;
        while ((pos = text.find(delimiter)) != std::string::npos) {
            tokens.push_back(std::stod(text.substr(0, pos)));
            text.erase(0, pos + delimiter.length());
        }
        tokens.push_back(std::stod(text));
        return tokens;
    }

    /**
     * Generates a random double from 0 to 1.
     * @return a random double
     */
    inline double nextDouble() {
        std::mt19937_64 rng(seeder());
        return probabilityDistributionDouble(rng);
    }

    /**
     * Generates a random bool.
     * @return random bool
     */
    inline bool nextBool() {
        std::mt19937 rng(seeder());
        return probabilityDistributionInt(rng) == 1;
    }

    /**
     * Generates a random int.
     * @param min bottom limit
     * @param max top limit
     * @return a random integer
     */
    inline unsigned int nextUnsignedInt(unsigned int min, unsigned int max) {
        std::uniform_int_distribution intDist{min, max};
        std::mt19937 rng(seeder());
        return intDist(rng);
    }
}


#endif //NEUROEVOLUTION_UTIL_H
