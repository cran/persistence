/**
 * @file eccezioni.h
 */

#pragma once

#include <stdexcept>
#include <string>
#include <format>
#include <source_location> // <-- Nuovo header C++20

namespace Milano {
    
    class eccezioni : public std::runtime_error {
    private:
        std::string msg;
        
    public:
        /**
         * @brief Costruisce l'eccezione catturando automaticamente file e riga.
         * @param arg Il messaggio d'errore.
         * @param loc Struttura C++20 che di default cattura la posizione della chiamata.
         */
        eccezioni(const std::string& arg,
                  const std::source_location& loc = std::source_location::current())
        : std::runtime_error(arg) {
            // Estrae file e riga dal parametro loc
            msg = std::format("{}:{}: {}", loc.file_name(), loc.line(), arg);
        }
        
        virtual ~eccezioni() noexcept = default;
        
        [[nodiscard]] const char* what() const noexcept override {
            return msg.c_str();
        }
    };
    
} // namespace Milano


