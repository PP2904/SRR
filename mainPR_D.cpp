#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>
#include <random>
#include <list>
#include <iomanip>
#include <stdlib.h>
#include "PR_D.h"
#include "rand_rounding.h"


//
// das Programm wird ausgeführt, falls restriction = 0
//

using namespace std;

//Proportional Response Dynamics

// returned:
// - Equilibrium prices
// - equilibrium allocations

//Paper: Distributed Algorithms via Gradient Descent for Fisher Markets


/*

int random_number(int lb, int ub) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 engine(seed);
    return (engine() % (ub - lb + 1)) + lb;
}
*/


//Main method
//PrDynamicsRestZero(num_bidders, num_goods, bidders_PRD, initPrices,
//                       num_iterations);
vector<double> PrDynamicsRestZero(int num_bidders_PRD, int num_goods, vector<Bidder> bidders_PRD, vector<double> initPrices,
               int num_iterations) {

    vector<double> utility(num_bidders_PRD);


    //FOR SCHLEIFE FÜR ANZAHL WIEDERHOLUNGEN DES GESAMTEXPERIMENTS
    //for (int iter = 0; iter < num_iterations; iter++) {

    //vector<Bidder> bidders_PRD(num_bidders_PRD);

    for (int k = 0; k < num_bidders_PRD; ++k) {
        bidders_PRD[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        for (auto &v: bidders_PRD[k].valuation) v = (random_number(1, 11) + random_number(1, 15));
        bidders_PRD[k].budget = random_number(1, 11) + random_number(1, 31);
        bidders_PRD[k].spent.resize(num_goods, bidders_PRD[0].budget / (double) num_goods);
    }


    //int num_iterations = 2000;
    vector<double> prices(num_goods);
    for (int it = 0; it < num_iterations; ++it) {

        //in jeder iteration werden die preise des guts i auf die menge der preise,
        // die jeder bidder ausgegeben hat, gesetzt
        for (int j = 0; j < num_goods; ++j) {
            prices[j] = 0;
            for (int i = 0; i < bidders_PRD.size(); ++i)
                prices[j] += bidders_PRD[i].spent[j];

        }
        //update der valuations und spents pro bidder
        vector<vector<double>> update(bidders_PRD.size(), vector<double>(num_goods)); //
        for (int i = 0; i < bidders_PRD.size(); ++i) {
            for (int j = 0; j < num_goods; ++j) {
                update[i][j] = bidders_PRD[i].valuation[j] * bidders_PRD[i].spent[j] / prices[j];

            }
        }

        //new bid vector for next iteration
        for (int i = 0; i < bidders_PRD.size(); ++i) {
            for (int j = 0; j < num_goods; ++j) {
                bidders_PRD[i].spent[j] =
                        bidders_PRD[i].budget * update[i][j] /
                        accumulate(update[i].begin(), update[i].end(), 0.0);

            }
        }

        //print für jeden bidder und jede iteration dessen allocation von Gut 0 bis n
        cout << "Iteration " << it << ":\n";
        for (int i = 0; i < bidders_PRD.size(); ++i) {
            cout << "Bidder " << i << ": " << bidders_PRD[i] << endl;
        }
        cout << endl;

    }

    //von Max utility und utility (im equilibrium sind diese gleich)


    for (int b = 0; b < num_bidders_PRD; ++b) {
        for (int i = 0; i < num_goods; ++i) {
            if (prices[i] == 0) {
                printf("prices is 0");
                exit(EXIT_FAILURE);
            }
            utility[b] += bidders_PRD[b].valuation[i] * bidders_PRD[b].spent[i] /
                          prices[i]; //Aufpassen wenn prices[i] = 0!
        }

    }


    cout << endl;
    cout << "Fraktionales/optimales Ergebnis: ";
    cout << endl;
    for (int i = 0; i < num_bidders_PRD; ++i) {
        cout << "Max Utility: " << std::setprecision(3) << utility[i] << endl;
    }
    cout << "\n";
    for (int i = 0; i < num_goods; ++i) {
        cout << prices[i] << "\n";
    }

    return initPrices;
}


