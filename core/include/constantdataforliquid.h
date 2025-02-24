// constants.h
#pragma once
#include <map>
#include <vector>
#include <string>

extern std::map<std::string, std::vector<double>> sigma_p_values_liquid_map;
extern std::map<std::string, std::vector<double>> sigma_p_values_liquidcoh_map;
//extern std::map<std::string, std::vector<std::vector<double>>> sigma_ex_values_map;



void initialize_atomic_dataforliquid() ; // 在头文件中声明初始化原子数据的函数
