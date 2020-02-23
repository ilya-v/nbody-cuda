#include "../params.hpp"

#include <cassert>
#include <iostream>


void test_params() {
	{
		params_t p = params_t::default_params();
		auto s = p.to_string();
		assert(std::string::npos != s.find("#n_steps              = 1000000"));
		assert(std::string::npos != s.find("#random_seed          = 102030"));
		assert(std::string::npos != s.find("#r2_eps               = 0.01"));
	}

	{
		params_t p;
		p.read("../../test/params-test-1.txt");
		auto s = p.to_string();
		assert(std::string::npos != s.find("#random_seed          = 2"));
		assert(std::string::npos != s.find("#config_r             = 23.5"));
		assert(std::string::npos != s.find("#dt_start             = 0.002"));
	}
}