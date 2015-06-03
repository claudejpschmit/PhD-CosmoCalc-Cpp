#pragma once

#include <map>
#include <string>
#include <vector>

using namespace std;

class CAMB_CALLER {

    public:
        CAMB_CALLER();
        ~CAMB_CALLER();
        void call(map<string, double> params);
    private:
        vector<string> file_content;
        string parameter_names [8];
        void update_params_ini(map<string, double> params);
        void create_output_file();
};
