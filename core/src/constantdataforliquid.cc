
#include "constantdataforliquid.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
// std::map<std::string, std::vector<double>> sigma_p_values_map;
// std::map<std::string, std::vector<double>> sigma_i_values_map;
// std::map<std::string, std::vector<std::vector<double>>> sigma_ex_values_map;

 std::map<std::string, std::vector<double>> sigma_p_values_liquid_map;
 std::map<std::string, std::vector<double>> sigma_p_values_liquidcoh_map;

void initialize_atomic_dataforliquid() {
    // 初始化不同原子的数据



 std::ifstream file("datawithcore/Lne-energytransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

    std::string line;
   std::vector<double> values_vector;
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquid_map["ne"] = values_vector; // 假设所有的数据都是CO2的
    }

  file.close();



file.open("datawithcore/Lne-momentumtransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

   values_vector.clear();
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquidcoh_map["ne"] = values_vector; // 假设所有的数据都是CO2的
    }
  file.close();

file.open("datawithcore/Lar-energytransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

   values_vector.clear();
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
       

        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquid_map["ar"] = values_vector; // 假设所有的数据都是CO2的
    }
 file.close();

file.open("datawithcore/Lar-momentumtransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

   values_vector.clear();
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
      

        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquidcoh_map["ar"] = values_vector; // 
    }
 file.close();

   // sigma_p_values_liquid_map["ar"] = {};  

   // sigma_p_values_liquidcoh_map["ar"] = {};  // 初始化对应的sigma_i_values
    

file.open("datawithcore/Lkr-energytransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

   values_vector.clear();
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
      

        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquid_map["kr"] = values_vector; // 假设所有的数据都是CO2的
    }
 file.close();

file.open("datawithcore/Lkr-momentumtransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

    values_vector.clear();
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
     

        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquidcoh_map["kr"] = values_vector; // 
    }
 file.close();



file.open("datawithcore/Lxe-energytransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

    values_vector.clear();
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
       

        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquid_map["xe"] = values_vector; // 假设所有的数据都是CO2的
    }
 file.close();

file.open("datawithcore/Lxe-momentumtransfer.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

   values_vector.clear();
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
 

        while (getline(iss, value, ',')) {
            values_vector.push_back(std::stod(value));
        }

        sigma_p_values_liquidcoh_map["xe"] = values_vector; // 
    }
 file.close();
//    sigma_p_values_liquid_map["kr"] = {
// };  // 初始化对应的sigma_p_values
//    sigma_p_values_liquidcoh_map["kr"] = {
// };  // 初始化对应的sigma_i_values
    

//    sigma_p_values_liquid_map["xe2"] = {
// };  // 初始化对应的sigma_p_values
//    sigma_p_values_liquidcoh_map["xe2"] = {
// };  // 初始化对应的sigma_i_values

   
   
}
