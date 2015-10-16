#include "Analyser.hpp"

Analyser::Analyser()
{}

Analyser::~Analyser()
{}

Fisher_return_pair Analyser::build_Fisher_inverse(vector<string> param_keys, string run_prefix,\
        string path)
{
    Fisher_return_pair RESULT;
    struct F_values 
    {
        string key1, key2;
        double value;
    };
    vector<F_values> F_ab;

    int num_params = param_keys.size();
    for (int i = 0; i < num_params; i++)
    {
        for (int j = i; j < num_params; j++)
        {

            //generate the relevant filenames
            string filename = path + run_prefix + "_Fisher_" + param_keys[i] +\
                              "_" + param_keys[j] + ".dat";
            ifstream f(filename);
            if (!f.is_open())
                filename = path + run_prefix + "_Fisher_" + param_keys[j] +\
                           "_" + param_keys[i] + ".dat";

            //First order file
            stringstream command_buff;
            command_buff << "python OrderFile.py " << filename;
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            system(command);
            delete command;

            //Read in the data
            ifstream file;
            file.open(filename);
            string line;
            vector<int> l;
            vector<double> F_l;
            while (getline(file,line))
            {
                int col1;
                double col2;
                istringstream ss(line);
                ss >> col1 >> col2;
                l.push_back(col1);
                F_l.push_back(col2);
            }
            file.close();

            //Then, construct the Fisher element F_key1_key2
            F_values F_ab_value;
            F_ab_value.key1 = param_keys[i];
            F_ab_value.key2 = param_keys[j];
            double v = 0;
            v = (2*l[0] + 1) * F_l[0];
            //This is probably quite inefficient, but it works
            /*
               for (int k = l[0]; k <= l[l.size()-1]; k++)
               {
               int index_low, index_high;
               for (int m = 0; m < l.size() - 1; m++)
               {
               if (k >= l[m] && k <= l[m+1])
               {
               index_low = m;
               index_high = m+1;
               break;
               }
               }
               double x1 = l[index_low];
               double x2 = l[index_high];
               double dx = x2-x1; 
               double y1 = F_l[index_low];
               double y2 = F_l[index_high];
               double dy = y2-y1;
            // interpolate linearly
            double fk = (dy/dx)*k + (x2*y1-x1*y2)/dx;  
            v += (2*k + 1) * fk;
            }
            */
            //here we sum up all F_l contributions - check what the right formula for this is.
            //int l_difference = l[1] - l[0];
            //for (int k = 0; k < l.size() - 1; k++)
            //{
            //    v += (2 * l[k] + 1) * F_l[k] * l_difference;
            //}
            //The last one we only add once.
            F_ab_value.value = v;
            F_ab.push_back(F_ab_value);
        }
    }
    //now we have all the necessary information in the F_ab vector
    //the only thing left is to put it in matrix form.
    vector<vector<vector<string>>> indecies;
    //size of the matrix is
    //int n = (-1+sqrt(1+8*filenames_Fl.size()))/2;
    mat F = randu<mat>(num_params,num_params);
    //fill the F matrix.
    for (int i = 0; i < num_params; i++)
    {
        vector<vector<string>> row;
        for (int j = 0; j < num_params; j++)
        {
            string key1, key2;
            vector<string> row_element;
            key1 = param_keys[i];
            key2 = param_keys[j];
            int F_ab_index = 0;
            for (int k = 0; k < F_ab.size(); k++)
            {
                if ((F_ab[k].key1 == key1 && F_ab[k].key2 == key2) ||\
                        (F_ab[k].key1 == key2 && F_ab[k].key2 == key1)){
                    F(i,j) = F_ab[k].value;
                    if (i<j){
                        row_element.push_back(key1);
                        row_element.push_back(key2);
                    }
                    else {
                        row_element.push_back(key2);
                        row_element.push_back(key1);
                    }

                }
            }

            row.push_back(row_element);
        }
        indecies.push_back(row);
    }
    /*
    cout << F << endl;
    vec eigenval;
    mat eigenvec;
    eig_sym(eigenval, eigenvec, F.i());
    cout << "e-values" << endl;
    for (int k = 0; k < 5; k++)
        cout << eigenval(k) << endl;
    */
    RESULT.matrix = pinv(F);
    RESULT.matrix_indecies = indecies;
    return RESULT;
}

Ellipse Analyser::find_error_ellipse(Fisher_return_pair finv, string param1, string param2, int run_number)
{
    int index1, index2;
    index1 = -1;
    index2 = -1;
    for (int i = 0; i < finv.matrix_indecies.size(); i++) {
        if ((index1 < 0) && (finv.matrix_indecies[0][i][1] == param1))  
            index1 = i;
        if ((index2 < 0) && (finv.matrix_indecies[0][i][1] == param2))  
            index2 = i;
    }
    if (index1 > index2){
        int buff = index1;
        index1 = index2;
        index2 = buff;
    }
    double sig_xx, sig_xy, sig_yy;
    sig_xx = finv.matrix(index1, index1);
    sig_xy = finv.matrix(index1, index2);
    sig_yy = finv.matrix(index2, index2);
    Ellipse ellipse;
    ellipse.a2 = (sig_xx + sig_yy)/2.0 + sqrt(pow(sig_xx - sig_yy,2)/4.0 + pow(sig_xy,2));
    ellipse.b2 = (sig_xx + sig_yy)/2.0 - sqrt(pow(sig_xx - sig_yy,2)/4.0 + pow(sig_xy,2));
    //theta is in radiants
    ellipse.theta = 0.5 * atan(2.0 * sig_xy/(sig_xx - sig_yy));
    ifstream runinfo("output/Fisher/RUN_INFO.dat");
    // find line ### run number bla ###
    //   while the line is not ### run number bla+1 ### or endoffile,
    //   check for the parameter keys and read the value into cx and cy.
    stringstream start_line, end_line;
    bool goodlines = false;
    start_line <<  "### run number " << run_number << " ###";
    end_line <<  "### run number " << run_number + 1 << " ###";
    while (runinfo.good())
    {   
        string line;
        getline(runinfo, line);
        if (line == start_line.str())
            goodlines = true;
        if (line == end_line.str())
            goodlines = false;
        if (goodlines)
        {
            if (line.find(param1) != string::npos)
            {
                string tmp1;
                char tmp2;
                double x;
                stringstream buff(line);
                buff >> tmp1 >> tmp2 >> x;
                ellipse.cx = x;
            }
            if (line.find(param2) != string::npos)
            {
                string tmp1;
                char tmp2;
                double y;
                stringstream buff(line);
                buff >> tmp1 >> tmp2 >> y;
                ellipse.cy = y;
            }

        }
    }

    runinfo.close();  
    return ellipse;
}

void Analyser::draw_error_ellipses(Fisher_return_pair finv, vector<string> param_keys, int run_number)
{
    // first need to know how many parameters we have, the grid size is equal to that
    int num_params = param_keys.size();
    // for each parameter pair we need to create an ellipse 
    vector<Ellipse> error_ellipses;
    for (int i = 0; i < num_params - 1; i++)
    {
        for (int j = i + 1; j < num_params; j++)
        {
            string param1 = param_keys[i];
            string param2 = param_keys[j];
            Ellipse ellipse = find_error_ellipse(finv, param2, param1, run_number);
            error_ellipses.push_back(ellipse);
        }
    }
    // a vector of ellipses should be passed to the drawing function
    //   as well as information about the corresponding parameters.
    // then draw
    string filename = "ellipse_info.tmp.dat";
    ofstream ellipse_file(filename);
    ellipse_file << num_params << endl;
    for (int i = 0; i<error_ellipses.size(); i++)
    {
        ellipse_file << error_ellipses[i].a2 << endl;
        ellipse_file << error_ellipses[i].b2 << endl;
        ellipse_file << error_ellipses[i].theta << endl;
        ellipse_file << error_ellipses[i].cx << endl;
        ellipse_file << error_ellipses[i].cy << endl;
    }
    
    //write a temporary file with the parameter names & pass to drawer
    ofstream param_file("paramfile.tmp.dat");
    for (int i = 0; i < num_params; i++)
    {
        param_file << param_keys[i] << endl;
    }
    param_file.close();
    stringstream command_buff;
    command_buff << "python plotEllipses.py " << filename << " paramfile.tmp.dat";
    char* command = new char[command_buff.str().length() + 1];
    strcpy(command, command_buff.str().c_str());
    system(command);
    delete command;
    system("rm ellipse_info.tmp.dat");
    system("rm paramfile.tmp.dat");   
}
