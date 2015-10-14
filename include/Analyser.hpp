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
        Fisher_return_pair build_Fisher_inverse(vector<string> filenames_Fl);
        Ellipse find_error_ellipse(Fisher_return_pair finv, string param1, string param2);
        void draw_error_ellipse(Ellipse ellipse, string param1, string param2);
};
