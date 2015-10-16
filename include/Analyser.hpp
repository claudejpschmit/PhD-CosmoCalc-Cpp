#include "Helper.hpp"
#include <armadillo>
#include <fstream>
#include <string>
#include "stdafx.h"
#include "interpolation.h"

using namespace alglib;
using namespace arma;

class Analyser {

    public:
        Analyser();
        ~Analyser();
        // TODO: These functions should be private and I
        // should have some public functions that do all of this at once.
        Fisher_return_pair build_Fisher_inverse(vector<string> param_keys, string run_prefix, string path);
        Ellipse find_error_ellipse(Fisher_return_pair finv, string param1, string param2, int run_number);
        void draw_error_ellipses(Fisher_return_pair finv, vector<string> param_keys, int run_number);
};
