#pragma once
#include "libomp/19.1.2/include/omp.h"
#include <array>
#include <vector>
#include <limits>
#include <complex>
#include <numeric>
#include <algorithm>
#include "Region.hpp"
#include "Rule.hpp"

/*
This behavior class is designed to integrate functions over a region utilizing a global adaptive strategy.
- 
*/

enum class IntegrationStrategy{
    parallel,
    serial,
};

template<int N>
class Integration {
    private:
        
        static const Rule<N>& m_get_rule(){
            static const Rule<N>& m_rule = Rule<N>::getInstance();
            return m_rule;
        }

        static const int& m_get_rule_func_calls(){
            static const int rule_func_call = m_get_rule().m_num_f_eval_total;
            return rule_func_call;
        }

        

        template<typename foo_return>
        static double error_estimation(std::array<foo_return, 6>& reduced_values){
            std::array<foo_return, 4> results_value;
            results_value.fill(foo_return(0.0));
            for (int i=0;i<6;i++){
                for (int j=0;j<4;j++){
                    results_value.at(j) += reduced_values.at(i) * m_get_rule().m_weights_rules.at( (j+1) ).at(i);
                }
            }
            double estimated_error{0.0};
            std::array<double, 3> null_list;
            for (int i=0; i<3; i++){
                null_list.at(i) = 0.0;
                for (int j=0; j<6; j++){
                    null_list.at(i) = std::max(null_list.at(i), std::abs(results_value.at(i+1)+ \
                       m_get_rule().m_scale_errors.at(i).at(j) * results_value.at(i)) * m_get_rule().m_norms.at(i).at(j) );
                }
            }
            if ( (m_get_rule().m_error_coefficiet.at(0) * null_list.at(0) <=  null_list.at(1)) &&\
                (m_get_rule().m_error_coefficiet.at(1) * null_list.at(1) <=  null_list.at(2)))
            {
                estimated_error = m_get_rule().m_error_coefficiet.at(2) * null_list.at(0);
            } else {
                double max_value = null_list.at(0);
                for (int i=1; i<3; i++){
                    if (max_value <= null_list.at(i)){
                        max_value = null_list.at(i);
                    }
                }
                estimated_error = m_get_rule().m_error_coefficiet.at(3) * max_value;
            }

            return estimated_error;
        }

        template<typename foo_t>
        static int get_direction( std::vector<foo_t>& evaluated_results, Region<N>& cur_region){
            std::array<double, N> fourth_difference_vec;
            // std::array<double, N> second_difference_vec;
            std::array<double, N> first_difference_vec;

            // std::vector<std::array<int, 3> > first_idx;
            // first_idx.push_back( std::array<int,3>{6*N+1, 6*N+4, 2*N+3} );
            // first_idx.push_back( std::array<int,3>{6*N+2, 6*N+3, 2*N+4} );
            // first_idx.push_back( std::array<int,3>{6*N+1, 6*N+2, 2*N+1} );
            // first_idx.push_back( std::array<int,3>{6*N+3, 6*N+4, 2*N+2} );

            for (int i=0; i<N; i++){
                fourth_difference_vec.at(i) = std::abs(evaluated_results.at(2*N+1+2*i)+evaluated_results.at(2*N+2+2*i)-2.0*evaluated_results.at(0)-\
                            m_get_rule().m_ratio_fourth_error *(evaluated_results.at(1+2*i)+evaluated_results.at(2+2*i)-2.0*evaluated_results.at(0)));
                if (fourth_difference_vec.at(i) <= 4.0 * std::numeric_limits<double>::epsilon() ) {
                    fourth_difference_vec.at(i) = 0.0;
                }

                // second_difference_vec.at(i) = std::abs( evaluated_results.at(6*N+1) + evaluated_results.at() - 2.0 * evaluated_results.at() \
                //                            + evaluated_results.at() + evaluated_results.at() - 2.0 * evaluated_results.at());
                // if (second_difference_vec.at(i) <= 4.0 * std::numeric_limits<double>::epsilon() ) {
                //     second_difference_vec.at(i) = 0.0;
                // }

                // if (N==2) {
                //     first_difference_vec.at(i) = std::abs(evaluated_results.at( first_idx.at(2*i).at(0) ) + evaluated_results.at(first_idx.at(2*i).at(1)) - 2.0 * evaluated_results.at(first_idx.at(2*i).at(2))) + \
                //                            std::abs( evaluated_results.at(first_idx.at(2*i+1).at(0)) + evaluated_results.at(first_idx.at(2*i+1).at(1)) - 2.0 * evaluated_results.at(first_idx.at(2*i+1).at(2)));
                //     if (first_difference_vec.at(i) <= 8.0 * std::numeric_limits<double>::epsilon() ) {
                //         first_difference_vec.at(i) = 0.0;
                //     }
                // }
                


            }
            double max_dir_value {0.0};
            int direction{0};
            for (int i=0; i<N; i++){
                if (fourth_difference_vec.at(i)>=max_dir_value){
                    max_dir_value = fourth_difference_vec.at(i);
                    direction = i;
                }
            }
            std::vector<int> index_set;
            for (int i=0; i<N; i++){
                if (max_dir_value == fourth_difference_vec.at(i)){
                    index_set.push_back(i);
                }
            }
            
            if (index_set.size()>1){
                double width_values_ref{ 0.0 };
                // You know in specific axis on which the function decay slowly, thus make the impact of the width less strong.
                // In general, when the decay rates of different axes are comparable, this is satisfactory for most applications.
                for (int i = 0; i<index_set.size(); i++){
                    double cur_scale = 1.0; // index_set.at(i)==1 ? 10.0 : 1.0;
                    if (cur_region.m_width.at(index_set.at(i))/cur_scale >= width_values_ref){
                        direction = index_set.at(i);
                        width_values_ref = cur_region.m_width.at(index_set.at(i));
                    }
                    
                }
            } else if (index_set.size() == 1) {
                direction = index_set.at(0);
            }
            return direction;
        }

        template <typename F>
        static auto call_rule(F f, Region<N>& cur_region) -> decltype(std::declval<F>()( std::declval< std::array<double,N> >() ) )
        {
            using foo_t = decltype( f( cur_region.m_center ) );
            std::vector<foo_t> evaluated_vec( m_get_rule_func_calls() );
            for (int i=0; i < m_get_rule_func_calls(); i++){
                std::array<double, N> cur_coordinate;
                for (int j=0; j<N; j++){
                    cur_coordinate.at(j) = cur_region.m_center.at(j) + 0.5 * m_get_rule().m_coordinates.at(i).at(j) * cur_region.m_width.at(j);
                }
                evaluated_vec.at(i) = f( cur_coordinate );
            }

            cur_region.divide_direction = get_direction(evaluated_vec, cur_region);
            cur_region.evaluated = true;

            std::array<foo_t, 6> reduced_values;
            reduced_values.fill(foo_t(0.0));
            reduced_values.at(0) = evaluated_vec.at(0);

            for(int i=0; i<2*N; i++){
                reduced_values.at(1) += evaluated_vec.at(1+i);
                reduced_values.at(2) += evaluated_vec.at(1+2*N+i);
                reduced_values.at(3) += evaluated_vec.at(1+4*N+i);
            }
            for (int i=0; i<(2*N*(N-1)); i++){reduced_values.at(4) += evaluated_vec.at(1+6*N+i);}
            for (int i=0; i<(1<<N); i++){reduced_values.at(5) += evaluated_vec.at(1+4*N+2*N*N+i);};

            double volumn{1.0};
            for (int i=0; i<N; i++){
                volumn *= cur_region.m_width.at(i)*0.5;
            }

            cur_region.estimated_error = error_estimation(reduced_values) * volumn;

            foo_t integral_val {0.0};
            for (int i=0; i<6; i++){
                integral_val += m_get_rule().m_weights_rules.at(0).at(i) * reduced_values.at(i);
            }
            if (std::abs(integral_val) < double(1<<N) * double(m_get_rule_func_calls()) * std::numeric_limits<double>::epsilon() )
            {
                integral_val = foo_t(0.0); 
            }
            integral_val = integral_val * volumn;
 
            return integral_val;
        }

    public:

    template<typename F>
    static auto integrate(F f, Region<N>& region, double& estimated_error, int& num_func_eval,  IntegrationStrategy integration_stategy=IntegrationStrategy::serial,\
                  double expected_error = 1e-6, int min_func_eval=500, int max_func_eval=20000)\
                ->decltype(std::declval<F>()( std::declval<std::array<double, N>>() ) )
    {

        num_func_eval = 0;
        using foo_t = decltype( f( region.m_center ) ); // It will be either real or complex double. 
        static_assert(!std::is_integral<foo_t>::value,
                "The return type must be either a real or complex floating point type.");
        using Region_Value = std::pair<Region<N> , foo_t>;
        std::vector<Region_Value> sub_regions_vector;

        foo_t estimated_integral = call_rule(f, region);
        sub_regions_vector.push_back( Region_Value(region, estimated_integral) );
        num_func_eval += m_get_rule_func_calls();

        while ( num_func_eval < max_func_eval )
        {
            switch (integration_stategy)
            {
                case IntegrationStrategy::parallel:
                    {
                        int number_current_threads = int(sub_regions_vector.size()) > m_get_rule().m_relaxation * omp_get_max_threads()? m_get_rule().m_relaxation * omp_get_max_threads() : int(sub_regions_vector.size()); 
                        #pragma omp parallel num_threads(number_current_threads)
                        {

                            Region_Value region_to_eval(region, estimated_integral);
                            #pragma omp critical
                            {
                                region_to_eval = sub_regions_vector.at(sub_regions_vector.size()-1);
                                sub_regions_vector.pop_back();
                            }
                            std::vector<Region<N> > splitted_sub_regions = region_to_eval.first.split();
                            std::array<foo_t, 2> subregion_integral;
                            #pragma omp barrier
                            for (int i=0; i<2; i++){
                                subregion_integral.at(i) = call_rule(f, splitted_sub_regions.at(i));
                            }
                            double error_global = std::abs(region_to_eval.second - subregion_integral.at(0) - subregion_integral.at(1));
                            if ((splitted_sub_regions.at(0).estimated_error + splitted_sub_regions.at(1).estimated_error) == 0){
                                // std::cout << "Here, left is " << splitted_sub_regions.at(0).estimated_error  << " right is " <<  splitted_sub_regions.at(1).estimated_error  << std::endl;
                                splitted_sub_regions.at(1).estimated_error = m_get_rule().m_error_coefficiet.at(5) * error_global;
                                splitted_sub_regions.at(0).estimated_error = m_get_rule().m_error_coefficiet.at(5) * error_global;

                            } else{
                                double error_tmp = m_get_rule().m_error_coefficiet.at(5)*error_global + splitted_sub_regions.at(0).estimated_error +\
                                            error_global * m_get_rule().m_error_coefficiet.at(4)*splitted_sub_regions.at(0).estimated_error\
                                            /(splitted_sub_regions.at(0).estimated_error + splitted_sub_regions.at(1).estimated_error);

                                splitted_sub_regions.at(1).estimated_error = m_get_rule().m_error_coefficiet.at(5)*error_global + splitted_sub_regions.at(1).estimated_error +\
                                                m_get_rule().m_error_coefficiet.at(4)*splitted_sub_regions.at(1).estimated_error * error_global\
                                                /(splitted_sub_regions.at(0).estimated_error + splitted_sub_regions.at(1).estimated_error);
                                splitted_sub_regions.at(0).estimated_error = error_tmp;
                            }
                            #pragma omp critical
                            {
                                for (int i=0; i<2; i++){
                                    sub_regions_vector.push_back( Region_Value(splitted_sub_regions.at(i), subregion_integral.at(i)));
                                }                
                            }
                        }
                        num_func_eval += 2 * m_get_rule_func_calls() * number_current_threads;
                        break;
                    }
                default:
                    {
                        Region_Value region_to_eval = sub_regions_vector.at(sub_regions_vector.size()-1);
                        sub_regions_vector.pop_back();
                        std::array<foo_t, 2> subregion_integral;
                        std::vector<Region<N> > splitted_sub_regions = region_to_eval.first.split();
                        for (int i=0; i<2; i++){
                            subregion_integral.at(i) = call_rule(f, splitted_sub_regions.at(i));
                        }
                        double error_global = std::abs(0.1*region_to_eval.first.estimated_error +  region_to_eval.second - subregion_integral.at(0) - subregion_integral.at(1));
                        if ((splitted_sub_regions.at(0).estimated_error + splitted_sub_regions.at(1).estimated_error) == 0){
                            // std::cout << "Here, left is " << splitted_sub_regions.at(0).estimated_error  << " right is " <<  splitted_sub_regions.at(1).estimated_error  << std::endl;
                            splitted_sub_regions.at(1).estimated_error = m_get_rule().m_error_coefficiet.at(5) * error_global;
                            splitted_sub_regions.at(0).estimated_error = m_get_rule().m_error_coefficiet.at(5) * error_global;

                        } else {
                            double error_tmp = m_get_rule().m_error_coefficiet.at(5)*error_global + splitted_sub_regions.at(0).estimated_error +\
                                        error_global * m_get_rule().m_error_coefficiet.at(4)*splitted_sub_regions.at(0).estimated_error\
                                        /(splitted_sub_regions.at(0).estimated_error + splitted_sub_regions.at(1).estimated_error);

                            splitted_sub_regions.at(1).estimated_error = m_get_rule().m_error_coefficiet.at(5)*error_global + splitted_sub_regions.at(1).estimated_error +\
                                            m_get_rule().m_error_coefficiet.at(4)*splitted_sub_regions.at(1).estimated_error * error_global\
                                            /(splitted_sub_regions.at(0).estimated_error + splitted_sub_regions.at(1).estimated_error);
                            splitted_sub_regions.at(0).estimated_error = error_tmp;
                        }


                        

                        for (int i=0; i<2; i++){
                            sub_regions_vector.push_back( Region_Value(splitted_sub_regions.at(i), subregion_integral.at(i)) );
                        }

                        num_func_eval += 2 * m_get_rule_func_calls();
                        break;
                }
            }

            estimated_integral = std::accumulate(sub_regions_vector.begin(), sub_regions_vector.end(), foo_t(0.0), \
                [](foo_t i, const Region_Value& cur_region_pair){
                    return i + cur_region_pair.second;});
            estimated_error = std::accumulate(sub_regions_vector.begin(), sub_regions_vector.end(), double(0.0), \
                [](double i, const Region_Value& cur_region_pair){
                    return i + std::abs(cur_region_pair.first.estimated_error);});


            if ( num_func_eval > min_func_eval ){
                if ( std::abs(estimated_error / estimated_integral) < expected_error)
                {
                    break;
                }
            }
            std::sort(sub_regions_vector.begin(), sub_regions_vector.end(), [](const Region_Value& lhs, const Region_Value& rhs){
                    return (std::abs(lhs.first.estimated_error) < std::abs(rhs.first.estimated_error)); });
        }
        return estimated_integral;
    }
};