#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <iomanip>
#include <stdlib.h>
#include <utility>
#include <algorithm>


using namespace std;

class Bidder {
public:
    vector<double> valuation; //was mir ein gut wert ist
    double budget;
    vector<double> spent; //für welches gut gibt bidder was aus (summe aller elem in spent ist budget)

    friend ostream &operator<<(ostream &os, const Bidder &b);
};

ostream &operator<<(ostream &os, const Bidder &b) {
    for (int j = 0; j < b.spent.size(); ++j) {
        os << b.spent[j] << " ";
    }
    return os;
}

int random_number(int lb, int ub) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 engine(seed);
    return (engine() % (ub - lb + 1)) + lb;
}


int SpendingRestrictedRound(vector<double> prices, vector<double> quantItemVec, double spendingRestriction,vector<Bidder> bidders, int num_goods, double valMultiplier, int num_bidders, int num_iterations, int num_iter_exp){

    //number of interested bidder in same good j (=interGood vector)
    vector<int> interGood(num_goods);

    //initialize interGood vector
    for (int j = 0; j < num_goods; ++j) {
        interGood[j] = 0;
    }


    //tuple: (number good, mbb)
    typedef pair<double, int> mytuple;
    double mbb = 0.0;
    //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...
    vector<vector<mytuple>> mbbVec(num_bidders, vector<mytuple>(num_goods));

    for (int i = 0; i < num_bidders; ++i) {
        cout << "Bidder: " << i << "\n";
        for (int j = 0; j < num_goods; ++j) {

            /* //test: 1 bidder bewertet das künstlich veränderte teuerste Gut als besonders wertvoll
             if (i == (num_bidders - (num_bidders - 1)) && j == (num_goods - 1)) {
                 mbb = 100;
                 mbbVec[i][j] = make_pair(mbb, j);
                 continue;
             }*/

            mbb = bidders[i].valuation[j] / prices[j];

            //interGood := Summe der Interessenten pro Gut
            if (mbb != 0.0) {
                interGood[j] += 1;
            }

            mbbVec[i][j] = make_pair(mbb, j);
            // cout << "Bidder " << i << " MBB for Good " << j << " is: " << bidders[i].valuation[j]/prices[j] << "\n";
        }
        //sortiere mbb vector nach größe bzgl. mbb wert (first value im tuple)
        sort(mbbVec[i].begin(), mbbVec[i].end(), greater<>());

        //for debugging
        //laufen über erste Reihe des mbbVec vectors und schauen alle Tuple an (p läuft durch Tuple in der ersten Reihe)
        for (const mytuple &p: mbbVec[i]) {
            cout << "Gut " << p.second << " hat MBB von: " << p.first << " und kostet " << prices[p.second] << "\n";
            //mbbAllocVec[i][p.second] = p.first;
        }


    }



    //for debugging
    cout << "\n";
    cout << "Die MBB Elemente in absteigender Reihenfolge pro Bidder: \n";
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            //const iterator (läuft nur über das aktuelle mbbVec[i] = aktuelle Reihe i)
            for (const mytuple &p: mbbVec[i]) {
                cout << "(Element " << p.second << " mit MBB " << p.first << " )" << " ";
            }
            //sehr wichtig!
            break;
        }
        cout << "\n";
    }

    cout << "\n";
    for (int i = 0; i < num_bidders; i++) {
        cout << "Für Bidder " << i << " greatest MBB " << "\n";
        for (int j = 0; j < num_goods; ++j) {

            const mytuple &p = mbbVec[i][j];

            cout << " for good " << p.second << " is: " << p.first << "\n";

        }

    }

    //debugging
    cout << "\n";
    for (int j = 0; j < num_goods; ++j) {
        cout << "Für Gut " << j << " gibt es " << interGood[j] << " Interessenten " << "\n";
    }
    cout << "\n";


    //maximum spend per item
    vector<double> spendPerItem(num_goods);
    for (int j = 0; j < num_goods; ++j) {
        spendPerItem[j] = 0;
    }


    //alloc mbb items
    vector<vector<double>> mbbItemAllocVec(num_bidders, vector<double>(num_goods));

    //spending vector  for spending graph Q(x)
    vector<vector<double>> spendVec(num_bidders, vector<double>(num_goods));

    //TODO: ich muss über mbbVec gehen und pro bidder das größte Gut zuweisen unter den Beschränkungen, die ich bereits verwendet habe
    // TODO: Ich benötige einen Durchlauf bis alle Items verkauft sind & jeder Bidder sein gesamtes Budget ausgegeben hat (!)





    //FIXME: Hier ist noch ein Fehler; die restriction wird overall nicht eingehalten (!)




        for (int iter = 0; iter < num_bidders; iter++) {

            for (const mytuple &p: mbbVec[iter]) {

                //number des aktuellen Goods
                int numGood = p.second;


                //Beschränkungen
                if (((spendPerItem[p.second] + ((min(bidders[iter].budget / prices[p.second], quantItemVec[p.second])))*prices[p.second]) <= spendingRestriction) && quantItemVec[p.second] != 0 && bidders[iter].budget != 0) {

                    //allocation; min <= ist überhaupt noch genug des Guts für den Bidder da?

                    mbbItemAllocVec[iter][p.second] = min(bidders[iter].budget / prices[p.second], quantItemVec[p.second]); // /interGood[j]);

                    //spending
                    spendVec[iter][p.second] = mbbItemAllocVec[iter][p.second] * prices[p.second];

                    //spending per Item maximum of 1
                    spendPerItem[p.second] = spendPerItem[p.second] + (mbbItemAllocVec[iter][p.second] * prices[p.second]);

                    //item wurde vekauft und muss daher dezimiert werden
                    quantItemVec[p.second] = quantItemVec[p.second] - mbbItemAllocVec[iter][p.second];

                    //stimm das mit dem budget abzug so?
                    // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                    bidders[iter].budget = bidders[iter].budget - spendVec[iter][p.second];

                    break;

                }

                if (((spendPerItem[p.second] + ((min(bidders[iter].budget / prices[p.second], quantItemVec[p.second])))*prices[p.second]) >=
                     spendingRestriction)) {


                    const mytuple &p = mbbVec[iter][numGood];

                    //allocation; min <= ist überhaupt noch genug des Guts für den Bidder da?

                    mbbItemAllocVec[iter][p.second] = min(bidders[iter].budget / prices[p.second], quantItemVec[p.second]); // /interGood[j]);
                    //spending
                    spendVec[iter][p.second] = mbbItemAllocVec[iter][p.second] * prices[p.second];
                    //spending per Item maximum of 1
                    spendPerItem[p.second] = spendPerItem[p.second] + (mbbItemAllocVec[iter][p.second] * prices[p.second]);
                    //item wurde vekauft und muss daher dezimiert werden
                    quantItemVec[p.second] = quantItemVec[p.second] - mbbItemAllocVec[iter][p.second];
                    //stimm das mit dem budget abzug so?
                    // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                    bidders[iter].budget = bidders[iter].budget - spendVec[iter][p.second];

                    break;

                }


            }

            //debugging - printing spending vector
            cout << "Spending graph: \n";
            for (int i = 0; i < num_bidders; ++i) {
                for (int j = 0; j < num_goods; ++j) {
                    cout << spendVec[i][j] << " ";
                }
                cout << "\n";
            }

        }






    for (int iter = 0; iter < num_bidders; iter++) {
        if (bidders[iter].budget <= 0.1) bidders[iter].budget = 0;
    }

    // Printing the vectors
    cout << "\n";
    //debugging - printing Alloc vector
    cout << "MBB alloc graph: \n";
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout << mbbItemAllocVec[i][j] << " ";
        }
        cout << "\n";
    }

    cout << "\n";

    //debugging - printing spending vector
    cout << "Spending graph: \n";
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout << spendVec[i][j] << " ";
        }
        cout << "\n";
    }


    //debugging - all goods allocated?
    cout << "\n";
    cout << "Restmenge jedes Guts: \n";
    for (int j = 0; j < num_goods; ++j) {
        cout << quantItemVec[j] << " ";
    }

    //debugging - max 1 dollar spend per good?
    cout << "\n";
    cout << "Ausgegeben pro Gut: \n";
    for (int j = 0; j < num_goods; ++j) {
        cout << spendPerItem[j] << " ";
    }


    //debugging - all budget spend?

    cout << "\n";
    cout << "Restbudget jedes Bidders: \n";
    for (int i = 0; i < num_bidders; ++i) {

        cout << bidders[i].budget << " ";
    }

    //debugging - utility
    double util = 0.0;

    cout << "\n";
    cout << "Utility: \n";
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            util += mbbItemAllocVec[i][j] * bidders[i].valuation[j];
        }
        cout << util << "\n";
    }


    return 0;

}


vector<double> ComputeFracEquilibrium(double spendingRestriction, vector<Bidder> bidders, int num_goods, double valMultiplier, int num_bidders, int num_iterations, double quantItem,
                            int num_iter_exp) {

    //prices vector
    vector<double> prices(num_goods);

    //FOR SCHLEIFE FÜR ANZAHL WIEDERHOLUNGEN DES GESAMTEXPERIMENTS
    for (int iter = 0; iter < num_iter_exp; iter++) {



        for (int it = 0; it < num_iterations; ++it) {

            //in jeder iteration werden die preise des guts i auf die menge der preise,
            // die jeder bidder ausgegeben hat, gesetzt
            for (int j = 0; j < num_goods; ++j) {
                prices[j] = 0;
                for (int i = 0; i < bidders.size(); ++i)
                    prices[j] += bidders[i].spent[j];

            }
            //update der valuations und spents pro bidder
            vector<vector<double>> update(bidders.size(), vector<double>(num_goods)); //
            for (int i = 0; i < bidders.size(); ++i) {
                for (int j = 0; j < num_goods; ++j) {
                    update[i][j] = bidders[i].valuation[j] * bidders[i].spent[j] / prices[j];

                }
            }

            //new bid vector for next iteration
            for (int i = 0; i < bidders.size(); ++i) {
                for (int j = 0; j < num_goods; ++j) {
                    bidders[i].spent[j] =
                            bidders[i].budget * update[i][j] / accumulate(update[i].begin(), update[i].end(), 0.0);

                }
            }

            //print für jeden bidder und jede iteration dessen allocation von Gut 0 bis n
            cout << "Iteration " << it << ":\n";
            for (int i = 0; i < bidders.size(); ++i) {
                cout << "Bidder " << i << ": " << bidders[i] << endl;
            }
            cout << endl;

        }

        //von Max utility und utility (im equilibrium sind diese gleich)

        vector<double> utility(num_bidders);
        vector<double> max_utility(num_bidders);
        for (int b = 0; b < num_bidders; ++b) {
            max_utility[b] = 0;
            for (int i = 0; i < num_goods; ++i) {
                if (prices[i] == 0) {
                    printf("prices is 0");
                    exit(EXIT_FAILURE);
                }
                utility[b] += bidders[b].valuation[i] * bidders[b].spent[i] / prices[i]; //Aufpassen wenn prices[i] = 0!
                if (max_utility[b] < bidders[b].valuation[i] / prices[i]) {
                    max_utility[b] = bidders[b].valuation[i] / prices[i];
                }
            }

            max_utility[b] *= bidders[b].budget;
        }

        // save utility from start
        vector<double> val_start(num_bidders);
        for (int b = 0; b < num_bidders; ++b) {
            for (int i = 0; i < num_goods; ++i) {
                val_start[b] = bidders[b].valuation[i];
            }
        }





        //set precision
        int pre = 3;

        //Optimales Ergebnis//

        cout << endl;
        cout << "Fraktionales/optimales Ergebnis: ";
        cout << endl;

        for (int i = 0; i < num_bidders; ++i) {
            cout << "Max Utility: " << std::setprecision(pre) << max_utility[i] << endl;
        }

        cout << endl;
        cout << "Bidders budget: ";
        cout << endl;

        for (int i = 0; i < num_bidders; ++i) {
            double rest = bidders[i].budget-accumulate(bidders[i].spent.begin(),bidders[i].spent.end(),0.0);
            cout << "restliches Budget bidder " << i  << ": " << std::setprecision(pre) << rest << endl;
        }

        cout << endl;

        //debugging:
        cout << endl;
        cout << "Allocation: ";
        cout << endl;

        vector<vector<double>> graph(num_bidders, vector<double>(num_goods));
        for (int i = 0; i < num_bidders; ++i) {
            for (int j = 0; j < num_goods; ++j) {
                graph[i][j] = bidders[i].spent[j] / prices[j]; //bidders.spent = nan
                if (isnan(graph[i][j])) {
                    graph[i][j] = 0.0000001;
                }
                cout << graph[i][j] << " ";
            }

            cout << "\n";

        }


       /* *//*** Write allocations to graph ***//*
        vector<vector<double>> graph(num_bidders, vector<double>(num_goods));
        for (int i = 0; i < num_bidders; ++i) {
            for (int j = 0; j < num_goods; ++j) {
                graph[i][j] = bidders[i].spent[j] / prices[j]; //bidders.spent = nan?
                if (isnan(graph[i][j])) {
                    graph[i][j] = 0.0000001;
                }


                //checken auf nan
                if (isnan(graph[i][j])) {
                    printf("graph value is nan");
                    exit(EXIT_FAILURE);
                }


            }
        }*/


    }

    return prices;
}


int main() {

    /********************
     *
     */

    //generate goods
    int num_goods;
    cout << "Number Goods: ";
    cin >> num_goods;

    //multiplier for valuation
    double valMultiplier = 1.;

    //generate bidders with val, budget and spent_vec randomly
    int num_bidders;
    cout << "Number Bidders: ";
    cin >> num_bidders;

    //laut paper muss number bidders < num goods sein
    if (num_bidders > num_goods) {
        printf("Error number bidders larger than number goods");
        exit(EXIT_FAILURE);
    }

    //spending restriction
    double spendingRestriction;
    cout << "What is the spending restriction?: ";
    cin >> spendingRestriction;


    //num_iterations = Anzahl der Iterationen des Handels auf dem FM
    int num_iterations;
    cout << "Number Iterations: ";
    cin >> num_iterations;

    //Quantität pro Gut
    double quantItem;
    cout << "Quantität eines Guts: ";
    cin >> quantItem;

    //quanitity of item i intially
    vector<double> quantItemVec(num_goods);
    for (int j = 0; j < num_goods; ++j) {
        quantItemVec[j] = quantItem;
    }

    //num_iter_exp = Anzahl Ausführungen des Gesamtexperiments
    int num_iter_exp;
    cout << "Number Iterations Experiment: ";
    cin >> num_iter_exp;

    vector<Bidder> bidders(num_bidders);

    for (int k = 0; k < num_bidders; ++k) {
        bidders[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        for (auto &v: bidders[k].valuation) v = (random_number(1, 11) + random_number(1, 15)) * valMultiplier;
        bidders[k].budget = random_number(1, 11) + random_number(1, 31);
        bidders[k].spent.resize(num_goods, bidders[0].budget / (double) num_goods);
    }


    //prices vector
    vector<double> prices(num_goods);

    /********************
     *
     */


    auto start = std::chrono::system_clock::now();


    //ofstream myfile;
    ofstream myfile2;


    myfile2.open("results.txt", std::ios_base::app);
    myfile2 << "Number Goods: " << num_goods << ", " << " Number Bidders: " << num_bidders << ", "
            << " Number Iterations: " << num_iterations << ", " << " Quantitaet pro Gut: " << quantItem << "\n";
    myfile2 << "max_utility for rounded alloc | max_utility" << "\n";

    /********************
    *                  *
    ********************/

    prices = ComputeFracEquilibrium(spendingRestriction,bidders, num_goods, valMultiplier, num_bidders, num_iterations, quantItem, num_iter_exp);

    SpendingRestrictedRound(prices, quantItemVec, spendingRestriction, bidders, num_goods, valMultiplier, num_bidders, num_iterations, num_iter_exp);

    /********************
    *                  *
    ********************/



    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    cout << "\n";
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;

}