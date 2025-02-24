/**
* @file noble gas/liquid
 * @author Cao Lei (caolei@ihep.ac.cn)
 * @version 1.0
 * @date 2025年12月23日
 *
 * @brief
 * Description:：This is a C++ program designed for electron transport in noble gases He,
 * Ne, Ar, Kr, or Xe, and liquids Ar, Kr, or Xe. The simulation, implemented in C++, is based on electron-atom collisions.
 *
 *
 * @section Features
 * - 1：gas He,Ne,Ar,Kr,Xe
 * - 2：liquid Ar,Kr,Xe
 *
 * @section description
 * 1. Compilation: Use CMake or compile directly via the command line with `g++`.
 * 2. Execution: Run the generated executable and follow the on-screen prompts.
 * 3. Configuration: Modify the  command-line parameters as needed.
 *
 *
 * @section Contact
 * If you have any questions or suggestions, please contact us via:
 * - Email:caolei@ihep.ac.cn
 *
 *
 * @section Acknowledgements
 * The gaseous cross-sections are from LXCat.
 */




//#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <cassert>
#include <numeric>
#include <cstdlib>
#include <cstddef>
#include "include/constantdata.h"
#include "include/constantdataforliquid.h"
#include <iomanip>
#include <omp.h>
#include <thread>
#include <atomic>
//#define M_PI 3.14159265

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//prepare for simulation

// Represents an electron with its properties and behavior
struct Electron {
    double x, y, z, t;
    std::vector<double> v; // Velocity
    std::vector<double> position; // Position coordinates
    double energy; // Kinetic energy of the electron
    double total_time; // Total time of simulation for this electron
    double time_step; // Time step for simulation
    double scatteringangle; // Scattering angle after collision
    bool uniform_field_ex; // Flag for uniform electric field
    double frequency; //
    double freemeantime; //
    double lambda; //
    double accel; // Acceleration due to electric field
    double dist; // Distance traveled

    // Constructor to initialize electron properties
    Electron(double total_time, double volts, std::vector<double> starting_pos,
             std::vector<double> v, std::mt19937 &generator, bool uniform_field_ex,
             double lambda,
             double frequency, double freemeantime) {
        // Initialize the electron
        double m_e = 9.1e-31; // Electron mass in kg
        double e = 1.60218e-19; // Elementary charge in C
        this->total_time = 0;
        this->v = v;
        this->energy = (0.5 * m_e * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])* 6.242e18);
        this->position = starting_pos;
        this->uniform_field_ex = uniform_field_ex;
        this->freemeantime = freemeantime;
        this->frequency = frequency;
        this->lambda = lambda;
        this->accel = volts * (e / m_e) * 100; // Calculate acceleration
    }
};


// Converts energy in electron volts to velocity
double E_to_v(double E_eV) {
    double v = std::sqrt((2 * E_eV * 1.60218e-19) / 9.1e-31); // m/s
    return v;
}

// Converts energy from Joules to electron volts
double J_to_eV(double J) {
    return J * 6.242e18;
}

// Calculates the Euclidean norm (magnitude) of a vector
double norm(std::vector<double> v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


// Generates a random unit vector
std::vector<double> random_unit_vector(std::mt19937 &gen) {
    std::normal_distribution<double> dist(0.0, 1.0);
    std::vector<double> v(3);
    for (int i = 0; i < 3; i++) {
        v[i] = dist(gen);
    }
    double norm_v = norm(v);
    if (norm_v != 0) {
        for (int i = 0; i < 3; i++) {
            v[i] /= norm_v;
        }
    }
    return v;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//for sampling
// Generates a random velocity vector based on a given temperature
std::vector<double> random_velocity(double T, double k_b, double M_a, std::mt19937 &gen) {
    double a = std::sqrt(k_b * T / M_a);
    std::normal_distribution<double> distribution(0.0, a);
    double vx = distribution(gen);
    double vy = distribution(gen);
    double vz = distribution(gen);
    return {vx, vy, vz};
}

int index(double v, double eV_min, double eV_max, int eV_number) {
    double m_e = 9.1e-31;
    double E = J_to_eV(0.5 * m_e * v * v);
    if (E < eV_min) {
        E = eV_min;
    }
    int ind = static_cast<int>(std::round((eV_number - 1) / log10(eV_max / eV_min) * log10(E / eV_min)));
    if (ind >= 1000) {
        ind = 999;
    }
    if (ind < 1) {
        ind = 1;
    }
    return ind;
}


// get cross-section
double x_sec_p(double v, const std::vector<double> &sigma, double eV_min, double eV_max, int eV_number) {
    int idx = index(v, eV_min, eV_max, eV_number);
    return sigma[idx];
}


// Retrieves the excitation cross-section for a given energy level and excitation level
double x_sec_ex(double v, double level, const std::vector<std::vector<double>> &sigma_ex, double eV_min, double eV_max, int eV_number) {
    int idx = index(v, eV_min, eV_max, eV_number);
    int level_int = static_cast<int>(level);
    double x_sec = 0.0;
    switch (level_int) {
        case 1:
            x_sec = sigma_ex[0][idx];
            break;
        case 2:
            x_sec = sigma_ex[1][idx];
            break;
        case 3:
            x_sec = sigma_ex[2][idx];
            break;
        default:
            x_sec = 0.0;
            break;
    }
    return x_sec;
}

double probability( double x_sec, double xsec_tot ,int factor_null_collision) {
    return  x_sec / xsec_tot / factor_null_collision;
}


// Simulates the time until the next collision
double get_meanfreetime(Electron &thiselec, std::mt19937 &gen) {//next_collision
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double random_number = dis(gen);
    double time_step = -log(random_number) / thiselec.frequency;
    if (time_step > 1.2 * thiselec.freemeantime) {
        time_step = 1.2 * thiselec.freemeantime;
    }
    if (time_step < 0.8 * thiselec.freemeantime) {
        time_step = 0.8 * thiselec.freemeantime;
    }

    if (time_step > 1e-11) {
        time_step = 1e-11;
    }
    // Cap the time step to prevent unreasonably large values
//    if (time_step > 1e-10) {
//        time_step = 1e-10;
//    }

    return time_step;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//particle change

// Updates the position and velocity of an electron
Electron update_pos_vel(Electron &Electron) {
    // Update position using the Verlet integration method
    Electron.position[2] += Electron.v[2] * Electron.time_step + 0.5 * Electron.accel * Electron.time_step * Electron.time_step;
    Electron.v[2] += Electron.accel * Electron.time_step;
    Electron.position[1] += Electron.v[1] * Electron.time_step;
    Electron.position[0] += Electron.v[0] * Electron.time_step;

    return Electron;
}


// Decreases the electron's energy and normalizes its velocity
void excitation(Electron &thiselec, double eV, double m_e) {
    thiselec.energy -= eV;
    double v_norm = sqrt((2 * thiselec.energy / 6.242e18) / m_e);
    for (int i = 0; i < 3; i++) {
        thiselec.v[i] = v_norm * (thiselec.v[i] / v_norm); // Normalize velocity
    }
}


// for gas. Simulates an elastic collision between an electron and another particle
Electron gas_elastic_collision(Electron &obj, double u, std::vector<double> vm, double m_e, double M_a, std::mt19937 gen) {
    std::vector<double> unit_vec = random_unit_vector(gen);
    std::vector<double> a(3);
    for (int i = 0; i < 3; i++) {
        a[i] = M_a * u * unit_vec[i];
    }
    std::vector<double> b(3);
    for (int i = 0; i < 3; i++) {
        b[i] = m_e * obj.v[i];
    }
    std::vector<double> c(3);
    for (int i = 0; i < 3; i++) {
        c[i] = M_a * vm[i];
    }
    std::vector<double> v1(3);
    for (int i = 0; i < 3; i++) {
        v1[i] = (a[i] + b[i] + c[i]) / (m_e + M_a);
    }

    double norm_v1 = 0.0;
    double norm_v2 = 0.0;
    double dot_product = 0.0;
    for (int i = 0; i < 3; i++) {
        norm_v1 += v1[i] * v1[i];
        norm_v2 += obj.v[i] * obj.v[i];
        dot_product += v1[i] * obj.v[i];
    }
    double angle_rad = acos(dot_product / (sqrt(norm_v1) * sqrt(norm_v2)));
    obj.scatteringangle = angle_rad * 180.0 / M_PI;
    obj.v = v1;
    return obj;
}


// Simulates ionization of an electron
void ionization(Electron &thiselec, int &total_ionizations, double volts,double ionization_energy,std::mt19937 generator) {

      //   std::vector v = random_unit_vector(generator) * E_to_v(1);
      // //  electron_list.push_back(new Electron(thiselec.t,  volts, thiselec.position, v,generator, thiselec.uniform_field_ex, thiselec.lambda, thiselec.frequency,thiselec.freemeantime));
      //   thiselec.energy=thiselec.energy-ionization_energy;
      //   total_ionizations++;


}


// Simulates an energy collision between an electron and another particle
void liquid_energy_collision(Electron &obj, double u, std::vector<double> vm, double m_e, double M_a, double col_prob_e, double col_prob_m, std::mt19937 gen) {
    // Store the original velocity for later use
    std::vector<double> origin_v = obj.v;

    // Generate a random unit vector for the direction of the collision
    std::vector<double> unit_vec = random_unit_vector(gen);

    // Calculate the momentum vectors for the electron and the other particle
    std::vector<double> a(3);
    for (int i = 0; i < 3; i++) {
        a[i] = M_a * u * unit_vec[i];
    }

    std::vector<double> b(3);
    for (int i = 0; i < 3; i++) {
        b[i] = m_e * obj.v[i];
    }

    std::vector<double> c(3);
    for (int i = 0; i < 3; i++) {
        c[i] = M_a * vm[i];
    }

    // Calculate the new velocity vector after the collision
    std::vector<double> v1(3);
    for (int i = 0; i < 3; i++) {
        v1[i] = (a[i] + b[i] + c[i]) / (m_e + M_a);
    }

    // Calculate the norms of the new and original velocity vectors
    double norm_v1 = 0.0;
    double norm_v2 = 0.0;
    double dot_product = 0.0;
    for (int i = 0; i < 3; i++) {
        norm_v1 += v1[i] * v1[i];
        norm_v2 += origin_v[i] * origin_v[i];
        dot_product += v1[i] * origin_v[i];
    }

    // Determine the new velocity based on the collision probabilities
    double r = col_prob_m / col_prob_e;
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    if (col_prob_m < col_prob_e) {
        if (dis(gen) > r) {
            // Scale the original velocity by the ratio of the norms
            for (int i = 0; i < 3; i++) {
                obj.v[i] = (sqrt(norm_v1) / sqrt(norm_v2)) * origin_v[i];
            }
        } else {
            obj.v = v1;
        }
    } else {
        if (dis(gen) < (1 / r)) {
            // Scale the new velocity by the ratio of the norms
            for (int i = 0; i < 3; i++) {
                obj.v[i] = (sqrt(norm_v2) / sqrt(norm_v1)) * v1[i];
            }
        } else {
            obj.v = v1;
        }
    }

    // Calculate the new scattering angle
    double norm_new = 0.0;
    for (int i = 0; i < 3; i++) {
        dot_product += origin_v[i] * obj.v[i];
        norm_new += obj.v[i] * obj.v[i];
    }

    if (norm_v2 == 0 || norm_new == 0) {
        std::cout << "Calculation is wrong" << std::endl;
    }

    double angle_rad = std::acos(dot_product / (sqrt(norm_new) * sqrt(norm_v2)));
    obj.scatteringangle = angle_rad * 180.0 / M_PI; // Convert radians to degrees
}









/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//update struction

void gas_update(Electron &thiselec, int &total_ionizations, double T,
                std::vector<double> sigma_p_values,
                std::vector<double> sigma_i_values, std::vector<std::vector<double>> sigma_ex,
                double eV_min, double eV_max, int eV_number, double m_e, double M_a, double k_b, double n,
                std::mt19937 &generator,int factor_null_collision,double ionization_energy,double volts) {
    // Calculate the time step until the next collision
    thiselec.time_step = get_meanfreetime(thiselec, generator);

    // Update the electron's position and velocity
    thiselec = update_pos_vel(thiselec);

    // Reset scattering angle and calculate the electron's energy
    thiselec.scatteringangle = 0;
    thiselec.energy = J_to_eV(0.5 * m_e * (thiselec.v[0] * thiselec.v[0] + thiselec.v[1] * thiselec.v[1] + thiselec.v[2] * thiselec.v[2]));

    // Generate a random velocity for the target particle
    std::vector<double> vm = random_velocity(T, k_b, M_a, generator);

    // Calculate the relative velocity magnitude between the electron and the target particle
    double u = norm({(thiselec.v[0] - vm[0]), (thiselec.v[1] - vm[1]), (thiselec.v[2] - vm[2])});

    // Retrieve cross-section values for various interactions

    double xsec_el = x_sec_p(norm(thiselec.v), sigma_p_values, eV_min, eV_max, eV_number);
    double xsec_ioni = x_sec_p(norm(thiselec.v), sigma_i_values, eV_min, eV_max, eV_number);
    double xsec_ex_1 = x_sec_ex(norm(thiselec.v), 11.55, sigma_ex, eV_min, eV_max, eV_number);
    double xsec_ex_2 = x_sec_ex(norm(thiselec.v), 13, sigma_ex, eV_min, eV_max, eV_number);
    double xsec_ex_3 = x_sec_ex(norm(thiselec.v), 14, sigma_ex, eV_min, eV_max, eV_number);
    double xsec_tot = xsec_el + xsec_ioni + xsec_ex_1 + xsec_ex_2 + xsec_ex_3;
//    std::cout<< "xsec_el "<< xsec_el <<std::endl;
//    std::cout<<"energy"<< thiselec.energy <<std::endl;
    // std::cout<< "xsec_tot "<< xsec_tot <<std::endl;
    //
//	thiselec.lambda=1/(xsec_tot * 1e-4)/(n * 1e6)/factor_null_collision;
//	thiselec.frequency=u/thiselec.lambda;
//	thiselec.freemeantime=1/thiselec.frequency;

    // Calculate probabilities for various interactions
    double col_prob = probability(xsec_el,xsec_tot,factor_null_collision);
    // std::cout<<"col_prob"<< col_prob <<std::endl;
    double ion_prob = probability( xsec_ioni,xsec_tot,factor_null_collision);
    double ex_1_prob = probability( xsec_ex_1,xsec_tot,factor_null_collision);
    double ex_2_prob = probability( xsec_ex_2,xsec_tot,factor_null_collision);
    double ex_3_prob = probability( xsec_ex_3,xsec_tot,factor_null_collision);

    // Randomly determine which interaction occurs
    std::uniform_real_distribution<double> rand_event(0.0, 1.0);
    double prob = rand_event(generator);
    bool hasinteracted = false;

    if (prob < ion_prob) {
       ionization(thiselec, total_ionizations,  volts, ionization_energy, generator);
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob)) {
        excitation(thiselec, 11.55, m_e);
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob + ex_2_prob)) {
        excitation(thiselec, 13, m_e);
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob + ex_2_prob + ex_3_prob)) {
        excitation(thiselec, 14, m_e);
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob + ex_2_prob + ex_3_prob + col_prob)) {
        gas_elastic_collision(thiselec, u, vm, m_e, M_a, generator);
        hasinteracted = true;
    }

    // Update the interaction distance and position


    // Update the electron's energy and total time
    thiselec.energy = J_to_eV(0.5 * m_e * (thiselec.v[0] * thiselec.v[0] + thiselec.v[1] * thiselec.v[1] + thiselec.v[2] * thiselec.v[2]));
    thiselec.total_time += thiselec.time_step;

    // Recalculate the target particle's velocity and relative velocity
    vm = random_velocity(T, k_b, M_a, generator);
    u = norm({(thiselec.v[0] - vm[0]), (thiselec.v[1] - vm[1]), (thiselec.v[2] - vm[2])});

    // Recalculate cross-section values for the updated conditions
    xsec_el = x_sec_p(norm(thiselec.v), sigma_p_values, eV_min, eV_max, eV_number);
    xsec_ioni = x_sec_p(norm(thiselec.v), sigma_i_values, eV_min, eV_max, eV_number);
    xsec_ex_1 = x_sec_ex(norm(thiselec.v), 11.55, sigma_ex, eV_min, eV_max, eV_number);
    xsec_ex_2 = x_sec_ex(norm(thiselec.v), 13, sigma_ex, eV_min, eV_max, eV_number);
    xsec_ex_3 = x_sec_ex(norm(thiselec.v), 14, sigma_ex, eV_min, eV_max, eV_number);
    xsec_tot = xsec_el + xsec_ioni + xsec_ex_1 + xsec_ex_2 + xsec_ex_3;

    // Update lambda for the new conditions
    thiselec.lambda=1/(xsec_tot * 1e-4 * n * 1e6 * factor_null_collision);
    thiselec.frequency=u/thiselec.lambda;
    thiselec.freemeantime=1/thiselec.frequency;
}









// Function to simulate updates to an electron's state based on various physical interactions
void liquid_update(Electron &thiselec, int &total_ionizations, double T,
                   std::vector<double> sigma_p_values_li, std::vector<double> sigma_p_values_liandcoh,
                   std::vector<double> sigma_i_values, std::vector<std::vector<double>> sigma_ex,
                   double eV_min, double eV_max, int eV_number, double m_e, double M_a, double k_b, double n,
                   std::mt19937 &generator,int factor_null_collision,double ionization_energy,double volts) {
    // Calculate the time step until the next collision
    thiselec.time_step = get_meanfreetime(thiselec, generator);

    // Update the electron's position and velocity
    thiselec = update_pos_vel(thiselec);

    // Reset scattering angle and calculate the electron's energy
    thiselec.scatteringangle = -1;
    thiselec.energy = J_to_eV(0.5 * m_e * (thiselec.v[0] * thiselec.v[0] + thiselec.v[1] * thiselec.v[1] + thiselec.v[2] * thiselec.v[2]));

    // Generate a random velocity for the target particle
    std::vector<double> vm = random_velocity(T, k_b, M_a, generator);

    // Calculate the relative velocity magnitude between the electron and the target particle
    double u = norm({(thiselec.v[0] - vm[0]), (thiselec.v[1] - vm[1]), (thiselec.v[2] - vm[2])});

    // Retrieve cross-section values for various interactions
    double xsec_el_e = x_sec_p(norm(thiselec.v), sigma_p_values_li, eV_min, eV_max, eV_number);

    double xsec_el_m = x_sec_p(norm(thiselec.v), sigma_p_values_liandcoh, eV_min, eV_max, eV_number);
    double xsec_ioni = x_sec_p(norm(thiselec.v), sigma_i_values, eV_min, eV_max, eV_number);
    double xsec_ex_1 = x_sec_ex(norm(thiselec.v), 11.55, sigma_ex, eV_min, eV_max, eV_number);

    double xsec_ex_2 = x_sec_ex(norm(thiselec.v), 13, sigma_ex, eV_min, eV_max, eV_number);
    double xsec_ex_3 = x_sec_ex(norm(thiselec.v), 14, sigma_ex, eV_min, eV_max, eV_number);
    double xsec_tot = xsec_el_e + xsec_el_m + xsec_ioni + xsec_ex_1 + xsec_ex_2 + xsec_ex_3;

    //
    thiselec.lambda=1/(xsec_tot * 1e-4)/(n * 1e6)/factor_null_collision;
    thiselec.frequency=u/thiselec.lambda;
    thiselec.freemeantime=1/thiselec.frequency;

    // Calculate probabilities for various interactions
    double col_prob_e = probability(xsec_el_e,xsec_tot,factor_null_collision);
    double col_prob_m = probability( xsec_el_m,xsec_tot,factor_null_collision);
    double ion_prob = probability( xsec_ioni,xsec_tot,factor_null_collision);
    double ex_1_prob = probability( xsec_ex_1,xsec_tot,factor_null_collision);
    double ex_2_prob = probability( xsec_ex_2,xsec_tot,factor_null_collision);
    double ex_3_prob = probability( xsec_ex_3,xsec_tot,factor_null_collision);

    // Randomly determine which interaction occurs
    std::uniform_real_distribution<double> rand_event(0.0, 1.0);
    double prob = rand_event(generator);
    bool hasinteracted = false;

    if (prob < ion_prob) {
        ionization(thiselec, total_ionizations,  volts, ionization_energy, generator);
        total_ionizations++;
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob)) {
        excitation(thiselec, 11.55, m_e);
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob + ex_2_prob)) {
        excitation(thiselec, 13, m_e);
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob + ex_2_prob + ex_3_prob)) {
        excitation(thiselec, 14, m_e);
        hasinteracted = true;
    } else if (prob < (ion_prob + ex_1_prob + ex_2_prob + ex_3_prob + col_prob_e + col_prob_m)) {
        liquid_energy_collision(thiselec, u, vm, m_e, M_a, col_prob_e, col_prob_m, generator);
        hasinteracted = true;
    }



    // Update the electron's energy and total time
    thiselec.energy = J_to_eV(0.5 * m_e * (thiselec.v[0] * thiselec.v[0] + thiselec.v[1] * thiselec.v[1] + thiselec.v[2] * thiselec.v[2]));
    thiselec.total_time += thiselec.time_step;

    // Recalculate the target particle's velocity and relative velocity
    vm = random_velocity(T, k_b, M_a, generator);
    u = norm({(thiselec.v[0] - vm[0]), (thiselec.v[1] - vm[1]), (thiselec.v[2] - vm[2])});

    // Recalculate cross-section values for the updated conditions
    xsec_el_e = x_sec_p(norm(thiselec.v), sigma_p_values_li, eV_min, eV_max, eV_number);
    xsec_el_m = x_sec_p(norm(thiselec.v), sigma_p_values_liandcoh, eV_min, eV_max, eV_number);
    xsec_ioni = x_sec_p(norm(thiselec.v), sigma_i_values, eV_min, eV_max, eV_number);
    xsec_ex_1 = x_sec_ex(norm(thiselec.v), 11.55, sigma_ex, eV_min, eV_max, eV_number);
    xsec_ex_2 = x_sec_ex(norm(thiselec.v), 13, sigma_ex, eV_min, eV_max, eV_number);
    xsec_ex_3 = x_sec_ex(norm(thiselec.v), 14, sigma_ex, eV_min, eV_max, eV_number);
    xsec_tot = xsec_el_e + xsec_el_m + xsec_ioni + xsec_ex_1 + xsec_ex_2 + xsec_ex_3;

    // Update lambda for the new conditions
    thiselec.lambda=1/(xsec_tot * 1e-4)/(n * 1e6)/factor_null_collision;
    thiselec.frequency=u/thiselec.lambda;
    thiselec.freemeantime=1/thiselec.frequency;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {//int argc, char** argv
    // 初始化MPI环境
//    MPI_Init(&argc, &argv);
//    int world_size;
//    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//
//    int world_rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
/*
    std::cout << "Please enter the electron energy in electron volts: eV: ";
    double elec_energy;
    std::cin >> elec_energy;

    std::cout << "Please enter the atom name(eg:ne, ar, kr, he, xe): ";
    std::string atom_name;
    std::cin >> atom_name;

    std::cout << "Please enter the electron emission position(x y z,separated by space): ";
    std::vector<double> pos(3);
    std::cin >> pos[0] >> pos[1] >> pos[2];
*/
    

    std::vector<int> voltages = {  100,  200,  500,  1000, 2000};
    std::vector<double> cutoffs= {  5E-4, 5E-4, 1E-3, 1E-3, 1E-3};
    int state=0;    // 如果是气态，为0，如果液态，为1
    int factor_null_collision=2;//gas:10,liquid:2
    std::string atom_name="ar";
    //if (world_rank == 0) {
    for (int v_index = 0; v_index < 1; v_index++)
    {
        int volts = voltages[v_index];
        double distance_cut = cutoffs[v_index];


//////////////////////////////////////////////////////
double elec_energy=0.1; //eV
    std::vector<double> pos(3);
    pos= {0,0,0};

    double M_a; // neon atom mass (kg)
    double T;

    std::vector<double> sigma_i_values;
    std::vector<std::vector<double>> sigma_ex;


    std::string defaultatom = "ar";
    initialize_atomic_data();
    sigma_i_values = sigma_i_values_map[defaultatom];
    sigma_ex = sigma_ex_values_map[defaultatom];
    std::vector<double> sigma_p_values_liandcoh;
    std::vector<double> sigma_p_values_li;
    std::vector<double> sigma_p_values;
    double ionization_energy;
    if(state==0){
        sigma_p_values= sigma_p_values_map[atom_name];
//                        for (double value : sigma_p_values) {
//                            std::cout << value << std::endl;
//                        }
    }

    if(state==1){
        initialize_atomic_dataforliquid();
        sigma_p_values_li = sigma_p_values_liquid_map[atom_name];
        sigma_p_values_liandcoh = sigma_p_values_liquidcoh_map[atom_name];
    }

    // 根据atom_name的值进行不同的赋值操作
    if (atom_name == "ar") {
        M_a = 6.63e-26;
        T = 87;
        ionization_energy=15.7;
    } else if (atom_name == "kr") {
        M_a = 13.6e-26;
        T = 120;
        ionization_energy=14;
    } else if (atom_name == "ne") {
        M_a = 3.35e-26;
        T = 27;
        ionization_energy=21.5;
    } else if (atom_name == "he") {
        M_a = 0.6646e-26;
        T = 4;
        ionization_energy=24.6;
    } else if (atom_name == "xe") {
        M_a = 21.802e-26;
        T = 165;
        ionization_energy=12.13;
    }
    double n;        // neon atomic density cm-3
    if(state==0){
        T=293;
        n=1e20;
    }else{ n = 1e22;}

    double m_e = 9.1e-31;   // electron mass (kg)
    double e = 1.60218e-19; // charge of the electron (in C)

    double eV_min = 0.0001; // (V)
    double eV_max = 200;
    int eV_number = 1000;
    double k_b = 1.38e-23; // Boltzmann constant (in J/K)
	////////////////////////////////////////////////////////////////////////////////////
     int threads = 100; // 
    omp_set_num_threads(threads);
   #pragma omp parallel for

        for (int i = 1; i < 100; i++){
			
			double lambda = 1e10, frequency = 1e11, freemeantime = 1e-11;
            bool uniform_field_ex = true;
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count() + rand()  );//omp_get_thread_num()
            std::vector<Electron> electron_list;
			std::vector<double> v = {random_unit_vector(generator)[0] * E_to_v(elec_energy), random_unit_vector(generator)[1] * E_to_v(elec_energy), random_unit_vector(generator)[2] * E_to_v(elec_energy)};
            electron_list.push_back(Electron(0, volts, pos, v, generator, uniform_field_ex,lambda,frequency, freemeantime));
			

                                         
                                         std::string folderPath = "result/";
                                         std::string filename;
                                         if(state==0){
                                             filename = folderPath + "gas-"+ std::to_string(factor_null_collision) + "-"+ atom_name +"_1e20_" + std::to_string(volts) + "V_" + std::to_string(i+1) + ".txt";
                                         }
                                         if(state==1){
                                             filename = folderPath + "liquid-"+ std::to_string(factor_null_collision) + "-"+ atom_name +"_1e20_" + std::to_string(volts) + "V_" + std::to_string(i+1) + ".txt";
                                         }

                                         std::ofstream fileID;
                                         fileID.open(filename, std::ios::app);

                                         double starting_z = 0;
    double t = electron_list[0].total_time = 0;
    double dt = electron_list[0].time_step;
    double x = electron_list[0].position[0];
    double y = electron_list[0].position[1];
    double z = electron_list[0].position[2];
    double energy = electron_list[0].energy;
    double vd = (z - starting_z) / t;
    double s = electron_list[0].scatteringangle;
    int total_ionizations = 0;
    int simulation_step = 1;
    const int write_every = 5000;
	while (z < distance_cut) { //electron_list[0].total_time < time_cut//electron_list[0].z<distance_cut
        Electron thiselec = electron_list[0];

        if (state==1){
            liquid_update(thiselec, total_ionizations, T, sigma_p_values_li, sigma_p_values_liandcoh, sigma_i_values, sigma_ex, eV_min, eV_max, eV_number, m_e, M_a, k_b, n, generator,factor_null_collision,ionization_energy,volts);

        }
        if (state==0){
            
            gas_update(thiselec, total_ionizations, T, sigma_p_values, sigma_i_values, sigma_ex, eV_min, eV_max, eV_number, m_e, M_a, k_b, n, generator,factor_null_collision,ionization_energy,volts);
          
        }
		 electron_list[0] = thiselec;
        simulation_step++;
        if (simulation_step % write_every != 0) {
            continue;
        }
            t = thiselec.total_time;
            dt = thiselec.time_step;
            double vx = thiselec.v[0];
            double vy = thiselec.v[1];
            double vz = thiselec.v[2];
            x = thiselec.position[0];
            y = thiselec.position[1];
            z = thiselec.position[2];
            energy = thiselec.energy;
            vd = std::abs((z - starting_z) / t);
			fileID << std::scientific;
            fileID << t * 1e9 << ", "  << x * 1e6 << ", " << y * 1e6 << ", " << z * 1e6 << ", "  << energy << ", " << vd  << "\n";
			
			}
      fileID.close();
                                     
     }

    }

    return 0;
}

