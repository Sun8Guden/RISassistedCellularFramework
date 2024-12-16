#pragma once
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <cassert>

/* 
This class is an exercise of using both "Strategy pattern" and "Singleton" pattern.

- Implementation state:
-   Only GenzMalik 7 is available. 
*/

enum class Rule_List {
    GenzMalik7
};


template <int N>
class Rule{
    double TWONDM { double(1<<N) }; 
    std::array<double, 6> m_lambdas;
    std::array<double, 4> m_lambdas_squares;

    void m_construct_lambda(){
        // Lambda is the generator related parameters.
        double delta_square = 0.47070000;
        double lbd1_square = 4.0 / ( 15.0 - 5.0/(delta_square) );
        double ratio_lbd1_delta = (1.0 - lbd1_square/delta_square) / 27.0;
        double lbd2_square = (5.0 - 7.0*lbd1_square - 35.0*ratio_lbd1_delta) / \
            (7.0 - 35.0*lbd1_square/3.0 - 35.0*ratio_lbd1_delta/delta_square);
        double lbdp_square = 0.5625;

        m_lambdas_squares = {lbd1_square, lbd2_square, lbdp_square, delta_square};
        m_lambdas = {0.0, std::sqrt(lbd2_square), std::sqrt(lbd1_square), std::sqrt(lbdp_square),\
                    std::sqrt(lbd1_square), std::sqrt(delta_square)};
        
        // for (int i=0; i<m_lambdas_squares.size(); i++){
        //     std::cout << m_lambdas_squares.at(i) << "  ";
        // }
        // std::cout << std::endl;
        // for (int i=0; i<m_lambdas.size(); i++){
        //     std::cout << m_lambdas.at(i) << "  ";
        // }
        // std::cout << std::endl;
        
    }

    void m_construct_coordinates(){
        std::array<double, N> origin_array;
        origin_array.fill(0.0);
        m_coordinates.push_back(origin_array);
        for (int i=0; i<N; i++){
            std::array<double, N> cur_array;
            cur_array.fill(0.0);
            cur_array.at(i) += m_lambdas.at(1);
            m_coordinates.push_back(cur_array);
            cur_array.at(i) -= 2.0 * m_lambdas.at(1);
            m_coordinates.push_back(cur_array);

        }
        for (int i=0; i<N; i++){
            std::array<double, N> cur_array;
            cur_array.fill(0.0);
            cur_array.at(i) += m_lambdas.at(2);
            m_coordinates.push_back(cur_array);
            cur_array.at(i) -= 2.0 * m_lambdas.at(2);
            m_coordinates.push_back(cur_array);
        }
        for (int i=0; i<N; i++){
            std::array<double, N> cur_array;
            cur_array.fill(0.0);
            cur_array.at(i) += m_lambdas.at(3);
            m_coordinates.push_back(cur_array);
            cur_array.at(i) -= 2.0 * m_lambdas.at(3);
            m_coordinates.push_back(cur_array);
        }
        for (int i=0; i<N-1; i++){
            for (int j=i+1; j<(N); j++){
                std::array<double, N> cur_array;
                cur_array.fill(0.0);
                cur_array.at(i) += m_lambdas.at(4);
                cur_array.at(j) += m_lambdas.at(4);
                m_coordinates.push_back(cur_array);
                cur_array.at(j) -= 2.0 * m_lambdas.at(4);
                m_coordinates.push_back(cur_array);
                cur_array.at(i) -= 2.0 * m_lambdas.at(4);
                m_coordinates.push_back(cur_array);
                cur_array.at(j) += 2.0 * m_lambdas.at(4);
                m_coordinates.push_back(cur_array);
            }
        }
        for (int i=0; i<(1<<N); i++){
            std::array<double, N> cur_array;
            cur_array.fill(0.0);
            int temp_i = i;
            for (int j=N-1; j>=0; j--){
                if (temp_i / (1<<j)){
                    cur_array.at(j) += m_lambdas.at(5);
                    temp_i -= (1<<j);
                } else {
                    cur_array.at(j) -= m_lambdas.at(5);
                }
            }
            m_coordinates.push_back(cur_array);
        }
        assert(m_coordinates.size() == m_num_f_eval_total );

        // for(int i=0; i<m_coordinates.size(); i++){
        //     for (int j=0; j<N; j++){
        //         std::cout << m_coordinates.at(i).at(j) << ", ";
        //     }
        //     std::cout<<std::endl;
        // }
    }

    void m_construct_weights_rules(){
        
        // Compute degree 7 rule weights. 
        m_weights_rules.at(0).at(5) = 1.0 / std::pow((3.0 * m_lambdas_squares.at(3)), 3.0) / TWONDM;
        m_weights_rules.at(0).at(4) = (1.0 - 5.0*m_lambdas_squares.at(3)/3.0) / (60.0*(m_lambdas_squares.at(0)-m_lambdas_squares.at(3))*m_lambdas_squares.at(0)*m_lambdas_squares.at(0));
        m_weights_rules.at(0).at(2) = (1.0 - 5.0*m_lambdas_squares.at(1)/3.0 - 5.0*TWONDM *m_weights_rules.at(0).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(1)))/\
                    (10.0 * m_lambdas_squares.at(0) * (m_lambdas_squares.at(0)-m_lambdas_squares.at(1))) - 2.0 * (double(N)-1.0) * m_weights_rules.at(0).at(4);
        m_weights_rules.at(0).at(1) = (1.0 - 5.0*m_lambdas_squares.at(0)/3.0 - 5.0*TWONDM*m_weights_rules.at(0).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(0)))/\
                    (10.0  *m_lambdas_squares.at(1) * (m_lambdas_squares.at(1)-m_lambdas_squares.at(0)));

        // Compute weight for degree 5 - 1; 
        m_weights_rules.at(1).at(5) = 1.0 / (36.0 * std::pow((m_lambdas_squares.at(3)), 3.0)) / TWONDM;
        m_weights_rules.at(1).at(4) = (1.0 - 9.0*TWONDM * m_weights_rules.at(1).at(5)*m_lambdas_squares.at(3)*m_lambdas_squares.at(3))/(36.0*m_lambdas_squares.at(0)*m_lambdas_squares.at(0));
        m_weights_rules.at(1).at(2) = (1.0 - 5.0*m_lambdas_squares.at(1)/3.0 - 5.0*TWONDM*m_weights_rules.at(1).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(1)))/\
                    (10.0 * m_lambdas_squares.at(0) * (m_lambdas_squares.at(0)-m_lambdas_squares.at(1))) - 2.0 * (double(N)-1.0) * m_weights_rules.at(1).at(4);
        m_weights_rules.at(1).at(1) = (1.0 - 5.0*m_lambdas_squares.at(0)/3.0 - 5.0*TWONDM*m_weights_rules.at(1).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(0)))/\
                    (10.0  *m_lambdas_squares.at(1) * (m_lambdas_squares.at(1)-m_lambdas_squares.at(0)));
    

        // Compute weight for degree 5 - 2; 
        m_weights_rules.at(2).at(5) = 5.0 / (108.0 * std::pow((m_lambdas_squares.at(3)), 3.0)) / TWONDM;
        m_weights_rules.at(2).at(4) = (1.0 - 9.0*TWONDM * m_weights_rules.at(2).at(5)*m_lambdas_squares.at(3)*m_lambdas_squares.at(3))/(36.0*m_lambdas_squares.at(0)*m_lambdas_squares.at(0));
        m_weights_rules.at(2).at(2) = (1.0 - 5.0*m_lambdas_squares.at(2)/3.0 - 5.0*TWONDM*m_weights_rules.at(2).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(2)))/\
                    (10.0 * m_lambdas_squares.at(0) * (m_lambdas_squares.at(0)-m_lambdas_squares.at(2))) - 2.0 * (double(N)-1.0) * m_weights_rules.at(2).at(4);
        m_weights_rules.at(2).at(3) = (1.0 - 5.0*m_lambdas_squares.at(0)/3.0 - 5.0*TWONDM*m_weights_rules.at(2).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(0)))/\
                    (10.0  *m_lambdas_squares.at(2) * (m_lambdas_squares.at(2)-m_lambdas_squares.at(0)));

        //Compute weight for degree 3
        m_weights_rules.at(3).at(5) = 1.0 / (54.0 * std::pow((m_lambdas_squares.at(3)), 3.0)) / TWONDM;
        m_weights_rules.at(3).at(4) = (1.0 - 18.0*TWONDM * m_weights_rules.at(3).at(5)*m_lambdas_squares.at(3)*m_lambdas_squares.at(3))/(72.0*m_lambdas_squares.at(0)*m_lambdas_squares.at(0));
        m_weights_rules.at(3).at(2) = (1.0 - 10.0*m_lambdas_squares.at(1)/3.0 - 10.0*TWONDM*m_weights_rules.at(3).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(1)))/\
                    (20.0 * m_lambdas_squares.at(0) * (m_lambdas_squares.at(0)-m_lambdas_squares.at(1))) - 2.0 * (double(N)-1.0) * m_weights_rules.at(3).at(4);
        m_weights_rules.at(3).at(1) = (1.0 - 10.0*m_lambdas_squares.at(0)/3.0 - 10.0*TWONDM*m_weights_rules.at(3).at(5) * m_lambdas_squares.at(3) * (m_lambdas_squares.at(3)-m_lambdas_squares.at(0)))/\
                    (20.0  *m_lambdas_squares.at(1) * (m_lambdas_squares.at(1)-m_lambdas_squares.at(0)));

        m_weights_rules.at(0).at(0) = TWONDM;


        for (int i=1; i<5; i++){
            for (int j=1; j<6; j++){
                m_weights_rules.at(i).at(j) = m_weights_rules.at(i).at(j) - m_weights_rules.at(0).at(j);
                m_weights_rules.at(i).at(0) = m_weights_rules.at(i).at(0) - double(m_num_f_eval_coordiate.at(j)) * m_weights_rules.at(i).at(j);
            }
        }

        for (int j=1; j<6; j++){
            m_weights_rules.at(0).at(j) = m_weights_rules.at(0).at(j) * TWONDM;
            m_weights_rules.at(0).at(0) = m_weights_rules.at(0).at(0) - double(m_num_f_eval_coordiate.at(j)) * m_weights_rules.at(0).at(j);
        }


    // for(int i=0; i<5; i++){
    //     for (int j=0; j<6; j++){
    //         std::cout << m_weights_rules.at(i).at(j) << ", ";
    //     }
    //     std::cout<<std::endl;
    // }


    }

    void m_construct_error_helpers(){
        // Scale and normalization for rules;
        m_ratio_fourth_error = m_lambdas.at(2)*m_lambdas.at(2)/m_lambdas.at(1)/m_lambdas.at(1);

        std::array<double, 6> WE;
        WE.fill(0.0);
        for (int i=0; i<3;i++){
            for (int j=0; j<6; j++){
                if (m_weights_rules.at(i+1).at(j)!=0.0)
                { m_scale_errors.at(i).at(j) = - m_weights_rules.at(i+2).at(j)/m_weights_rules.at(i+1).at(j);}
                else
                { m_scale_errors.at(i).at(j) = 100.0;}

                for (int k=0; k<6; k++){
                    WE.at(k) = m_weights_rules.at(i+2).at(k) + m_scale_errors.at(i).at(j) * m_weights_rules.at(i+1).at(k);
                }
                m_norms.at(i).at(j) = 0.0;
                for (int k=0; k<6; k++){
                    m_norms.at(i).at(j) += double(m_num_f_eval_coordiate.at(k)) * std::abs(WE.at(k));
                }
                m_norms.at(i).at(j) = double(1<<N) / m_norms.at(i).at(j);
            }
        }
        // for(int i=0; i<3; i++){
        //     for (int j=0; j<6; j++){
        //         std::cout << m_scale_errors.at(i).at(j) << ", ";
        //     }
        //     std::cout<<std::endl;
        // }
        // std::cout << "Section 2:\n";
        
        // for (int i=0; i<6; i++){
        //     std::cout << WE.at(i) << ", ";
        // }
        // std::cout<<std::endl;
        // for(int i=0; i<3; i++){
        //     for (int j=0; j<6; j++){
        //         std::cout << m_norms.at(i).at(j) << ", ";
        //     }
        //     std::cout<<std::endl;
        // }
    }

protected:
    // Construct the rule;
    Rule():m_weights_rules(std::vector<std::vector<double> >(5, std::vector<double>(6, 0.0))), \
            m_scale_errors(std::vector<std::vector<double> >(3, std::vector<double>(6, 0.0))), \
            m_norms(std::vector<std::vector<double> >(3, std::vector<double>(6, 0.0))) {
        m_construct_lambda();
        m_construct_coordinates();
        m_construct_weights_rules();
        m_construct_error_helpers();
        // std::cout << "The rule for " << N << " dimension integration is successfully initializated." << std::endl;
    }

public:

    std::vector<std::array<double, N> > m_coordinates;
    std::array<int, 6> m_num_f_eval_coordiate {1, 2*N, 2*N, 2*N, 2*N*(N-1), 1<<N};
    int m_num_f_eval_total { 1 + 4 * N + 2 * N * N + (1 << N)};
    std::vector<std::vector<double> > m_weights_rules;
    std::array<double, 6> m_error_coefficiet {5.0, 5.0, 1.0, 5.0, 0.5, 0.25};
    double m_ratio_fourth_error {0.0};
    std::vector<std::vector<double> > m_scale_errors;
    std::vector<std::vector<double> > m_norms;
    int m_relaxation { 3 };

    Rule(Rule& other_rule) = delete;
    void operator=(const Rule& other_rule) = delete;

    static Rule& getInstance(){
        static Rule rule;
        return rule;
    }
};