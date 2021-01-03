
//
// Created by Peter Pfrommer on 03.01.21.
//

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <chrono>
#include <random>
#include <list>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include "rand_rounding.h"
//#include "PR_D.h"


//
//
//this method rounds our SRR-Equilibrium allocation
//
//


using namespace std;


/*
int random_number(int lb, int ub) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 engine(seed);
    return (engine() % (ub - lb + 1)) + lb;
}

*/


//Main method
vector<vector<double>> roundingSRE(int num_bidders, int num_goods, vector<vector<double>> allocVec,
        int quant, vector<Bidder> bidders, vector<double> utilityRound ) {

    int pre = 3;

    //für budget printout (debugging)
    ofstream myfile2;
    myfile2.open("roundedResult.txt", std::ios_base::app);


    //fraktionale allokation von Bieter i und gut j
    vector<vector<double>> fractional_allocations(num_bidders, vector<double>(num_goods, 0.0));
    vector<vector<double>> integral_allocations(num_bidders, vector<double>(num_goods, 0.0));
    vector<vector<double>> final_allocations(num_bidders, vector<double>(num_goods, 0.0));


    //summe fraktionale Teile pro Gut über alle bidder
    vector<double> sum_frac(num_goods, 0.0);

    //Bieter 1 von Gut 1: 0,3 Bieter 2 v G1: 0,2 Bieter 3 v G1: 0,5
    //dann partial_sums = {0.3, 0.5, 1.0}
    //partial sums pro Gut (über alle Bidder)
    vector<double> partial_sums(num_bidders, 0.0);


    //die fraktionalen allokationen aus allocVec[i][j] werden auf fractional_allocations[i][j] addiert
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            fractional_allocations[i][j] = ((quant * (allocVec[i][j])) - (floor(quant * (allocVec[i][j]))));
            if (fractional_allocations[i][j] < 0.01 || fractional_allocations[i][j] > 0.99) {
                fractional_allocations[i][j] = 0;
            }
            //fractional_allocation[i][j]+integral_allocation[i][j] = allocVec[i][j]
            integral_allocations[i][j] = round(quant * allocVec[i][j] - fractional_allocations[i][j]);
            final_allocations[i][j] = integral_allocations[i][j];
            //cout << "Bidder " << i << " has " << fractional_allocations[i][j] << " of good " << j << "\n";
        }
    }

// pro bidder berechne ich die summe der fraktionalen Teile pro (jeweils ein) gut j; wird für partial sums benötigt;
//sum_frac ist damit die summe der fraktionalen Teile (über alle Bidder) eines guts j
    for (int j = 0; j < num_goods; ++j) {
        for (int i = 0; i < num_bidders; ++i) {
            //jeder fraktionale Teil eines jeden bidders i des spezifischen gutes j wird aufaddiert:
            sum_frac[j] = sum_frac[j] + fractional_allocations[i][j];

        }
        //cout << "Gut " << j << " hat in Summe " << sum_frac[j] << " fraktionale Einheiten" << "\n";
    }


    //zuweisung der fraktionalen Teile auf den (per rndm number) gezogenen Bidder
    for (int j = 0; j < num_goods; ++j) {
        //sum_frac[j] = summe fraktionale teile Gut j. sum_frac[j]=0 => kein fraktionaler teil
        //Wird gut fraktional aufgeteilt? Wenn nicht -> continue
        if (sum_frac[j] < 0.01) continue;
        //intitialisierung der partial sums mit den fraktionalen werten (nur so ist rekursion in for schleife über partial sums möglich
        partial_sums[0] = fractional_allocations[0][j];

        //startet bei i=1 (wegen rekursion)
        for (int i = 1; i < num_bidders; ++i) {
            partial_sums[i] = partial_sums[i - 1] + fractional_allocations[i][j];
        }
        for (int i = 0; i < num_bidders; ++i) {
            //sum_frac[j] is die summe der fraktionalen Allokationen pro gut j über alle Bidder i
            partial_sums[i] /= sum_frac[j];
        }

        //zufallszahl zw. 0 und 1 (double)
        //srand(time(NULL)); funktioniert nicht so gut

        double rdm_number = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        for (int i = 0; i < num_bidders; ++i) {
            //wenn zufallszahl <= partial_sums[i] => bieter i bekommt das fraktionale gut zugewiesen und break;
            if (rdm_number <= partial_sums[i]) {
                final_allocations[i][j] += sum_frac[j];
                break;
            }
        }
    }

    cout << "Original allocations:" << endl;
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            if ((quant * allocVec[i][j]) < 0.01) {
                allocVec[i][j] = 0;
            }
            cout << quant * allocVec[i][j] << " ";
        }
        cout << "|";
    }
    cout << endl;

    cout << "Fractional allocations:" << endl;
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout << fractional_allocations[i][j] << " ";
        }
        cout << "|";
    }
    cout << endl;

    cout << "Integral allocations:" << endl;
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout << integral_allocations[i][j] << " ";
        }
        cout << "|";
    }
    cout << endl;

    cout << "Randomized rounding allocations: " << endl;
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout << final_allocations[i][j] << " ";
        }
        cout << "|";
    }
    cout << endl;



    //die neuen/upgedateten/update max_utils berechnen und ausdrucken

    //utilityRound der gerundeten Alloks berechnen
    cout << "\n";
    cout << "utilityRound for rounded alloc | utilityRound: | integrality gap: \n";
    //myfile << "utilityRound for rounded alloc | utilityRound: | integrality gap: \n";
    //myfile2 << ", \n";


    double rd_util = 0.0;
    vector<double> rd_utilityRound(num_bidders);

    double int_gap = 0.0;
    double print_int_gap = 0.0;
    double avg_int_gap = 0.0;


    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            rd_util = rd_util + (((final_allocations[i][j]) / quant) * bidders[i].valuation[j]);
        }
        //utilityRound for rounded alloc
        rd_utilityRound[i] = rd_util;
        cout << rd_util << " | ";
        //myfile << rd_util << " | ";
        myfile2 << rd_util << " | ";

        //utilityRound:
        cout << std::setprecision(pre) << utilityRound[i] << " | ";
        //myfile << std::setprecision(pre)  << utilityRound[i] << " | ";
        myfile2 << std::setprecision(pre) << utilityRound[i] << "\n";

        //integrality gap:
        if (rd_utilityRound[i] <= utilityRound[i]) {
            cout << std::setprecision(pre) << rd_utilityRound[i] / utilityRound[i] << "\n";
            // myfile << std::setprecision(pre)  << rd_utilityRound[i]/utilityRound[i] << "\n";
            int_gap = int_gap + (rd_utilityRound[i] / utilityRound[i]);
        }
            //integer sol > optimal sol
        else {
            cout << " \n";
        }



    }
    return final_allocations;
}








