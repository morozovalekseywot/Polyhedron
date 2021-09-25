//
// Created by 159-egi on 2/16/17.
//

#ifndef INFINIMESH_CONFIGURATION_CONTROLLER_H
#define INFINIMESH_CONFIGURATION_CONTROLLER_H

#include <allstd.h>
#include <json11.hpp>

/**
 * @brief Представляет обертку для JSON файла
 */
class Configuration {
public:
    Configuration() = default;

    /**
     * @brief Базовый конструктор менеджера конфигурации
     * @param filename Путь к файлу конфигурации
     */
    explicit Configuration(string filename);

    /**
     * @brief Проверяет, существует ли в конфигурации параметр с данным именем
     * @param name Имя параметра для проверки на наличие в конфигурации
     * @return true если параметр существует в конфигурации, false в противном
     * случае
     */
    bool exist(const string &name) const;

    /**
     * @brief Проверяет, существует ли в конфигурации параметр с данным именем
     * в данной секции
     * @param section Раздел параметра для проверки на наличие в конфигурации
     * @param name Имя параметра для проверки на наличие в конфигурации
     * @return true если параметр существует в конфигурации, false в противном
     * случае
     */
    bool exist(const string &section,
               const string &name) const;

    /**
     * @brief Проверяет, существует ли в конфигурации параметр с данным именем
     * в данной секции
     * @param section Раздел параметра для проверки на наличие в конфигурации
     * @param subsection Подраздел раздела параметра для проверки на наличие в
     * конфигурации
     * @param name Имя параметра для проверки на наличие в конфигурации
     * @return true если параметр существует в конфигурации, false в противном
     * случае
     */
    bool exist(const string &section,
               const string &subsection,
               const string &name) const;


    /**
     * @brief Выбирает параметр из конфигурации
     * @param name Имя параметра, который выбирается из конфигурации
     * @return Возвращает obj[name], если параметр существует, иначе кидает
     * исключение
     */
    Configuration operator()(const string &name) const;

    /**
     * @details Выбирает параметр из конфигурации
     * @param section Раздел в котором содержится параметр
     * @param name Имя параметра
     * @return Возвращает obj[section][name], если параметр в секции существует,
     * иначе кидает исключение.
     */
    Configuration operator()(const string &section,
                             const string &name) const;

    /**
     * @details Выбирает параметр из конфигурации
     * @param section Раздел в котором содержится параметр
     * @param subsection Подраздел, в котором содержится параметр
     * @param name Имя параметра
     * @return Возвращает obj[section][subsection][name], если параметр
     * существует в выбранной секции и подсекции, иначе кидает исключение.
     */
    Configuration operator()(const string &section,
                             const string &subsection,
                             const string &name) const;

    /** Проверка на массив */
    bool is_array() const;

    /** Приводит JSON к типу boolean или кидает исключение */
    bool to_bool() const;

    /** Приводит JSON к типу integer или кидает исключение */
    int to_int() const;

    /** Приводит JSON к типу unsigned или кидает исключение */
    uint to_uint() const;

    /** Приводит JSON к типу double или кидает исключение */
    double to_double() const;

    /** Приводит JSON к типу string или кидает исключение */
    string to_string() const;

    /** Приводит JSON к массиву или кидает исключение */
    vector<Configuration> to_array() const;

    /// @brief Возвращает ассоциативный контейнер или кидает исключение
    std::map<string, double> to_map() const;

protected:
    /**
     * @brief Конструктор для создания временных объектов
     */
    explicit Configuration(json11::Json config_part);


public:
    /** Разобранный файл конфигурации */
    json11::Json m_config;
};

#endif //INFINIMESH_CONFIGURATION_CONTROLLER_H