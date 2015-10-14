#include "Analyser.hpp"

Analyser::Analyser()
{}

Analyser::~Analyser()
{}

Fisher_return_pair Analyser::build_Fisher_inverse(vector<string> filenames_Fl)
{
    Fisher_return_pair RESULT;
    struct F_values 
    {
        string key1, key2;
        double value;
    };
    //keep a list of parameter keys that have been used.
    vector<string> parameter_keys;
    vector<F_values> F_ab;
    for (int i = 0; i < filenames_Fl.size(); i++)
    {
        cout << filenames_Fl[i] << endl;
        //First order file
        stringstream command_buff;
        command_buff << "python OrderFile.py " << filenames_Fl[i];
        char* command = new char[command_buff.str().length() + 1];
        strcpy(command, command_buff.str().c_str());
        system(command);
        delete command;
        //Analyse filename so that we know what is in the file.
        int uscore_count = 0;
        string key1, key2;
        for (char & c : filenames_Fl[i])
        {
            if (c == '_' || c == '.')
                uscore_count++;
            if (uscore_count == 2 && c != '_')
                key1 += c;
            if (uscore_count == 3 && c != '_')
                key2 += c;
        }
        //Read in the data
        ifstream file;
        file.open(filenames_Fl[i]);
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
            cout << col1 << " " <<col2 << endl;
        }
        file.close();
        
        //Then, construct the Fisher element F_key1_key2
        F_values F_ab_value;
        bool k1new = true;
        bool k2new = true;
        for (int j = 0; j < parameter_keys.size(); j++)
        {
            if (parameter_keys[j] == key1)
                k1new = false;
            if (parameter_keys[j] == key2)
                k2new = false;
        }
        if (k1new)
            parameter_keys.push_back(key1);
        if (k2new)
            parameter_keys.push_back(key2);


        F_ab_value.key1 = key1;
        F_ab_value.key2 = key2;
        double v = 0;
        int l_difference = l[1] - l[0];
         //here we sum up all F_l contributions - check what the right formula for this is.
        for (int j = 0; j < l.size() - 1; j++)
        {
            v += (2 * l[j] + 1) * F_l[j] * l_difference;
        }
         //The last one we only add once.
        v += (2 * l[l.size()-1] + 1) * F_l[l.size()-1];
        F_ab_value.value = v;
        cout << v << endl; 
        F_ab.push_back(F_ab_value);
    }
    //now we have all the necessary information in the F_ab vector
    //the only thing left is to put it in matrix form.
    vector<vector<vector<string>>> indecies;
    //size of the matrix is
    //int n = (-1+sqrt(1+8*filenames_Fl.size()))/2;
    int n = parameter_keys.size();
    cout << n << endl;
    mat F = randu<mat>(n,n);
    //fill the F matrix.
    for (int i = 0; i < parameter_keys.size(); i++)
    {
        vector<vector<string>> row;
        for (int j = 0; j < parameter_keys.size();j++)
        {
            string key1, key2;
            vector<string> row_element;
            key1 = parameter_keys[i];
            key2 = parameter_keys[j];
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
    cout << "here" << endl;
    cout << F << endl;
    RESULT.matrix = pinv(F);
    RESULT.matrix_indecies = indecies;
    return RESULT;
}

Ellipse Analyser::find_error_ellipse(Fisher_return_pair finv, string param1, string param2)
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
    ellipse.theta = 0.5 * atan(2.0 * sig_xy/(sig_xx - sig_yy));
    return ellipse;
}

void Analyser::draw_error_ellipse(Ellipse ellipse, string param1, string param2)
{

}
