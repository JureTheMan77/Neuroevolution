//
// Created by jure on 7/15/20.
//

#ifndef NEUROEVOLUTION_LOGGING_H
#define NEUROEVOLUTION_LOGGING_H

#include <string>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace logging {

    /**
    * Gets current readable time. HH:MM:SS:millis
    */
    static std::string getCurrentTime() {
        std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
        std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

        std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
        std::time_t t = s.count();
        std::size_t fractional_seconds = ms.count() % 1000;

        std::ostringstream stream;
        stream << std::put_time(std::localtime(&t), "%H:%M:%S:") << fractional_seconds;

        return stream.str();
    }

    /**
     * Prints a @param message with a timestamp to stdout.
     * @param message message to print
     */
    static void logs(const std::string &message) {
        if (message.find('\n') != std::string::npos) {
            std::cout << "[" << getCurrentTime() << "] " << std::endl << message << std::endl;
        } else {
            std::cout << "[" << getCurrentTime() << "] " << message << std::endl;
        }
    }

    /**
     * Prints a @param message with a timestamp to stdout.
     * @param message message to print
     */
    [[maybe_unused]] static void logw(std::string message);

    /**
     * Prints a @param message with a timestamp to stderr.
     * @param message message to print
     */
    [[maybe_unused]] static void loge(std::string message);

}


#endif //NEUROEVOLUTION_LOGGING_H
