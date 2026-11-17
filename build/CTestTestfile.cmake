# CMake generated Testfile for 
# Source directory: /Users/guodongsun/Documents/GitHub/RISassistedCellularFramework
# Build directory: /Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test("option_pricing" "/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/build/examples/OptionPricingExample")
set_tests_properties("option_pricing" PROPERTIES  LABELS "validation" _BACKTRACE_TRIPLES "/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/CMakeLists.txt;67;add_test;/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/CMakeLists.txt;0;")
add_test("regression" "/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/build/tests/regression_test")
set_tests_properties("regression" PROPERTIES  _BACKTRACE_TRIPLES "/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/CMakeLists.txt;74;add_test;/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/CMakeLists.txt;0;")
add_test("monte_carlo_validation" "/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/build/validation/monte_carlo" "400000")
set_tests_properties("monte_carlo_validation" PROPERTIES  LABELS "validation" _BACKTRACE_TRIPLES "/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/CMakeLists.txt;84;add_test;/Users/guodongsun/Documents/GitHub/RISassistedCellularFramework/CMakeLists.txt;0;")
