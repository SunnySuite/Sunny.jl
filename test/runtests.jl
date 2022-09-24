""" To run tests, execute `test` in package mode with the Sunny package activated,
or, in the REPL, execute: `using Pkg; Pkg.test("Sunny")`. To execute only a single
test from the test suite, execute: `'Pkg.test("Sunny", test_args=["test_symmetry"])`,
for example, replacing `test_symmetry` with the name of the desired test.
"""

using TestItemRunner

@run_package_tests

