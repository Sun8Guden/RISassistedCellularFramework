#include <iostream> 
#include <functional> 
#include <memory>
#include <string>
#include <cxxabi.h>

template <class T>
std::string
type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
           (
#ifndef _MSC_VER
                abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                           nullptr, nullptr),
#else
                nullptr,
#endif
                std::free
           );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}

class FunctionTest{
    public:
    std::function<int(int)> m_fn;

    // FunctionTest(const std::function<int(int)> func):m_fn(func){
    //     std::cout << "non-reference Initialized; " << std::endl;
    // };

    // FunctionTest(const std::function<int(int)>& func):m_fn(func){
    //     std::cout << "reference Initialized; " << std::endl;
    // };

    FunctionTest(std::function<int(int)>&& func):m_fn(func){
        std::cout << "Move Initialized; " << std::endl;
    };

    // FunctionTest(std::function<int(int)>&& func):m_fn(func){
    //     std::cout << "Non-const, Move Initialized; " << std::endl;
    // };
};
 
void executeFunction(const std::function<int(int)>& func, int value) { 
    std::cout << "Again" << std::endl;
    func(value); 
} 
 
int main() { 
    double a {1.0};
    auto const myLambda = [&](int x) { std::cout << "Lambda called with value: " << x + a << std::endl; return x; }; 
     
    // executeFunction(myLambda, 5); // Output: Lambda called with value: 5 
     
    std::cout << type_name<decltype(myLambda)>() << std::endl;
    FunctionTest func_test(myLambda);
    func_test.m_fn(5);
    executeFunction(myLambda, 5);
    return 0; 
} 