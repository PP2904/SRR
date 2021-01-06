
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



/*int random_number(int lb, int ub) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 engine(seed);
    return (engine() % (ub - lb + 1)) + lb;
}*/




//Main method
vector<vector<double>> roundingSRE(int num_bidders, int num_goods, vector<vector<double>> &allocVec,
        int quant, vector<Bidder> &bidders, vector<double> MaxUtility, double spendingRestriction,int num_iterations, vector <double> &newPrices) {


    int pre = 3;

    //für budget printout (debugging)
    ofstream myfile2;
    myfile2.open("roundedResult.txt", std::ios_base::app);

    myfile2 << "Bidders: " << num_bidders << " Goods: " << num_goods << " Iterations: " << num_iterations
            << " spending restriction: " << spendingRestriction << " quantity per item: " << quant << "\n";


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
            fractional_allocations[i][j] = (((allocVec[i][j])) - (floor((allocVec[i][j]))));
            if (fractional_allocations[i][j] < 0.01 || fractional_allocations[i][j] > 0.99) {
                fractional_allocations[i][j] = 0;
            }
            //fractional_allocation[i][j]+integral_allocation[i][j] = allocVec[i][j]
            integral_allocations[i][j] = round(allocVec[i][j] - fractional_allocations[i][j]);
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

        //double rdm_number = (random_number(0, 5))/6;

        double rdm_number = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);


        for (int i = 0; i < num_bidders; ++i) {
            //wenn zufallszahl <= partial_sums[i] => bieter i bekommt das fraktionale gut zugewiesen und break;
              //if (rdm_number <= partial_sums[i]) {
                if (rdm_number <= partial_sums[i] && (bidders[i].budget - sum_frac[j]/newPrices[j]) <= 0){
                final_allocations[i][j] += sum_frac[j];
                break;
            }
        }
    }

    //debugging purpose
    //cout << "quantity per item is: " << quant << "\n";

    cout << "Original allocations:" << endl;
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            if ((allocVec[i][j]) < 0.01) {
                allocVec[i][j] = 0;
            }
            cout << allocVec[i][j] << " ";
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
           //cout << floor(final_allocations[i][j]) << " ";

           //round half up
           if(fractional_allocations[i][j] > 0.4){
               final_allocations[i][j] = ceil(fractional_allocations[i][j]);
               cout << ceil(final_allocations[i][j]) << " ";
           }
           else{
               final_allocations[i][j] = floor(fractional_allocations[i][j]);
               cout << floor(final_allocations[i][j]) << " ";
           }

        }
        cout << "|";
    }
    cout << endl;



    //die neuen/upgedateten/update max_utils berechnen und ausdrucken

    //utilityRound der gerundeten Alloks berechnen
    cout << "\n";
    cout << "utility for rounded alloc | utility: \n";
    myfile2 << "utility for rounded alloc | utility: \n";



    double rd_util = 0.0;
    vector<double> rd_utility(num_bidders);



    for (int i = 0; i < num_bidders; ++i) {
        rd_util = 0.0;
        for (int j = 0; j < num_goods; ++j) {
            if(final_allocations[i][j] == 1){
                rd_util += (((fractional_allocations[i][j])) * bidders[i].valuation[j]);
            }
            else{
                rd_util += 0.0;
            }
        }
        //utility for rounded alloc
        rd_utility[i] = rd_util;
        cout << rd_utility[i] << " | ";
        myfile2 << rd_utility[i] << " | ";

        //max utility
        cout << std::setprecision(pre) << MaxUtility[i] << " | " << "\n";
        myfile2 << std::setprecision(pre) << MaxUtility[i] << "\n";


    }
    return final_allocations;
}








