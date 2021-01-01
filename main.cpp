//
// Created by Peter Pfrommer on 19.12.20.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <list>
#include <iomanip>
#include <stdlib.h>
#include <utility>
#include <map>
#include <algorithm>
#include "mainPR_D.cpp"
#include "PR_D.h"


using namespace std;

/*
class Bidder {
public:
    vector<double> valuation; //was mir ein gut wert ist
    double budget;
    vector<double> spent; //für welches gut gibt bidder was aus (summe aller elem in spent ist budget)

    friend ostream &operator<<(ostream &os, const Bidder &b);
};
*/

ostream &operator<<(ostream &os, const Bidder &b) {
    for (int j = 0; j < b.spent.size(); ++j) {
        os << b.spent[j] << " ";
    }
    return os;
}

//type def für: vector<vector<mytuple>> mbbVec
typedef pair<double, int> myTuple;

/*

int random_number(int lb, int ub) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 engine(seed);
    return (engine() % (ub - lb + 1)) + lb;
}
*/


//return mbbVec; vector<vector<myTuple>> mbbVec(num_bidders, vector<myTuple>(num_goods));
//bspw: bidder 1: (mbb, 0), (mbb,1), ...
vector<vector<myTuple>> mbbGraph(int num_bidders, int num_goods, vector<Bidder> &bidders, vector<double> &newPrices,
                                 vector<vector<myTuple>> &mbbVec) {


    //sorted vector for mbb
    vector<vector<myTuple>> SortedMbbVec(num_bidders, vector<myTuple>(num_goods));

    double mbb;


    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            //FIXME: macht die if-Abfrage Sinn?
            //if(mbbVec[i][j].first > 0) {

            mbb = bidders[i].valuation[j] / newPrices[j];
            if (mbb < 0.001) mbb = 0;

            //Paar-Reihenfolge ist verkehrt zum Print-out
            mbbVec[i][j] = make_pair(mbb, j);
            //}
            //continue;

        }


        //sortiere mbb vector nach größe bzgl. mbb wert (first value im tuple)
        sort(mbbVec[i].begin(), mbbVec[i].end(), greater<>());
        SortedMbbVec = mbbVec;
    }

    mbbVec = SortedMbbVec;

    return SortedMbbVec;

}


// geht über den sortedVec und entscheidet pro Bidder, welche Güter er kauft;
// dabei darf maximal "spendingRestriction" GE auf jedes Gut verbraucht werden (in summe über alle Bidder)
vector<vector<double>> spendingGraph(int num_bidders, int num_goods, vector<Bidder> &bidders, vector<double> &newPrices,
                                     double spendingRestriction, vector<double> &spendPerItem,
                                     vector<double> &quantItem,
                                     vector<vector<myTuple>> &SortedMbbVec, vector<vector<double>> spendVec) {


    //spendVec = spending vector  for spending graph Q(x);
    // Bidder xy: Gut1, Gut2, ...



    //FIXME: spendPerItem wird bei jeder iteration iter wieder auf 0 gesetzt, warum?
    // => weil keine Referenz (&) übergeben wurde von spendPerItem ...



    //spending graph erzeugen
    for (int iter = 0; iter < num_bidders; iter++) {
        //laufen durch den sortierten Vektor
        //spending-graph edges are a subset of the mbb-graph edges
        //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...

        for (const myTuple &p: SortedMbbVec[iter]) {

            //number des aktuellen Goods
            int numGood = p.second;



            //new share of good
            double newShare = double(double(bidders[iter].budget / num_goods) / newPrices[p.second]);

            if (bidders[iter].budget == 0) {
                break;
            }


            if ((double(spendPerItem[p.second] + (newShare * newPrices[p.second])) <= spendingRestriction) &&
                quantItem[p.second] - newShare >= 0.0 && bidders[iter].budget != 0.0) {


                //spendVec wird erhöht durch neues share des Guts * price des guts
                spendVec[iter][p.second] += newShare * newPrices[p.second];


                //item wurde vekauft und muss daher dezimiert werden
                quantItem[p.second] -= newShare;
                if (quantItem[p.second] < 0.001) quantItem[p.second] = 0.0;

                //bisher für gut ausgegeben (overall agents)
                spendPerItem[p.second] += newShare * newPrices[p.second];


                // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                bidders[iter].budget -= newShare * newPrices[p.second];

                break;
            }

        }

        //spendPerItem = spendPerItem_Clone;

    }


    return spendVec;
}

vector<int> interestsGood(int num_bidders, int num_goods, vector<vector<myTuple>> &SortedMbbVec) {

    //wird jeden Schleifendurchlauf NEU berechnen

    //number of interested bidder in same good j (=interGood vector)
    vector<int> interGood(num_goods);


    for (int i = 0; i < num_bidders; ++i) {
        for (const myTuple &p: SortedMbbVec[i]) {

            if (p.first != 0.0) {
                interGood[p.second] += 1;
            }
        }

    }


    return interGood;
}


vector<double> currentPrice(int num_bidders, int num_goods, vector<Bidder> &bidders, vector<double> initPrices,
                            double spendingRestriction, vector<double> &spendPerItem, vector<double> &quantItem,
                            vector<vector<myTuple>> &SortedMbbVec,
                            vector<int> &interGood, vector<vector<double>> &spendVec, vector<vector<double>> &update,  vector<double> initBudget) {

    ofstream myfile;
    myfile.open("data.txt", std::ios_base::app);

    //passen Preise schrittweise an
    vector<double> newPrices(num_goods);


    newPrices = initPrices;



    //PR-D preisanpassung; warum sinken die mbbs dennoch auf 0 mit iterations -> unendlich (in dem Fall schon 25 =) ??
    // Vermutung: die Preise werden immer geringer und dadurch das Verhältnis von valuation/preis ... )

    for (int k = 0; k < num_goods; ++k) {
        //TODO: verändertert, dass alle Budgets ausgegeben werden (oder alle items aufgebraucht)
        newPrices[k] = 0;
        for (int i = 0; i < bidders.size(); ++i){
            newPrices[k] += bidders[i].spent[k];
        }
    }


    /* //idea: if interGood > bidders/2 => preise erhöhen
     //idea: if interGood <= bidders/2 => preise erniedrigen
     if(interGood[k] > (num_bidders/2)){
         newPrices[k] += bidders[i].spent[k];
     }
     else{
         newPrices[k] -= 100*(bidders[i].spent[k]);
     }*/



    //Problem ist, dass spent bereits im 2. Durchlauf = 0 ist. <- wird nicht upgedatet => update = 0 und dann teilen wir im
    //3. for loop hier durch update = 0 ...

    //erklärung: budget des bidders geht auf 0 (da er alles ausgegeben hat) und dann wird spent 0 => update = 0 und teilen durch 0 :(

    for (int i = 0; i < bidders.size(); ++i) {
        for (int j = 0; j < num_goods; ++j) {
            update[i][j] = (bidders[i].valuation[j] * bidders[i].spent[j]) / newPrices[j];

        }
    }


    //TODO: hier wird eingestiegen, wenn update vector 0 (fast 0) ist (!)
    for (int i = 0; i < bidders.size(); ++i) {
        for (int j = 0; j < num_goods; ++j) {
            //FIXME: Problem: irgendwann ist der update vector 0, weil spent-vector sehr klein ist und daher
            // nur noch wenig hinzukommt;
            if(accumulate(update[i].begin(), update[i].end(), 0.0) <= 0.001){
                cout << endl;
                //wieviel bleibt pro Gut übrig:
                cout << "available items: \n";
                for (int j = 0; j < num_goods; ++j) {
                    if (quantItem[j] < 0.01) {
                        quantItem[j] = 0;
                    }
                    cout << "Good " << j << " : " << quantItem[j] << "\n";
                    myfile << "Good " << j << " : " << quantItem[j] << "\n";
                }
                for (int i = 0; i < bidders.size(); ++i) {
                    cout << "(Bidder " << i << " (spend: " << 100*(1-(bidders[i].budget/initBudget[i])) << " %)"  << "\n";
                }
                //for debugging
                for (int j = 0; j < num_goods; ++j) {
                    cout << "\n";
                    cout << "Für Good " << j << " wurden " <<  spendPerItem[j] << " Geldeinheiten ausgegeben" << "\n";
                 }
                printf("update vector is zero");
                exit(EXIT_FAILURE);
            }
            bidders[i].spent[j] = ( bidders[i].budget * update[i][j] )  / accumulate(update[i].begin(), update[i].end(), 0.0);

        }
    }




/* FUNKTIONSAUFRUFE */

    SortedMbbVec = mbbGraph(num_bidders, num_goods, bidders, newPrices, SortedMbbVec);

    interGood = interestsGood(num_bidders, num_goods, SortedMbbVec);

    //TODO update spendVec bei jedem Loop
    spendVec = spendingGraph(num_bidders, num_goods, bidders, newPrices, spendingRestriction, spendPerItem, quantItem,
                             SortedMbbVec,
                             spendVec);


    //for debugging
    cout << "\n";
    cout << "mbb-graph: \n";
    for (int i = 0; i < num_bidders; ++i) {
        cout << "Bidder " << i << ": " << "\n";
        //const iterator (läuft nur über das aktuelle mbbVec[i] = aktuelle Reihe i)
        for (const myTuple &p: SortedMbbVec[i]) {
            //for debugging
            /*if (i == 0) {
                myfile << "(" << p.second << "," << setprecision(3) << p.first << ")" << " ";
            }*/
            cout << "(" << p.second << "," << setprecision(3) << p.first << ")" << " ";
        }
        cout << "\n";
        //myfile << "\n";
    }

    //for debugging
    cout << "\n";
    cout << "InterGood graph: \n";
    for (int j = 0; j < num_goods; ++j) {
        cout << interGood[j] << " ";
    }
    cout << "\n";


    /* //for debugging
     cout << "\n";
     cout << "Spending graph: \n";
     for (int i = 0; i < num_bidders; ++i) {
         for (int j = 0; j < num_goods; ++j) {
             cout << spendVec[i][j] << " ";
         }
         cout << "\n";
     }
 */

    return newPrices;

}

//Gleichgewichtspreise nach PR-Dynamics Algo => Vergleichswert
vector<double> PrDynamics(int num_bidders, int num_goods, vector<Bidder> &bidders_PRD, vector<double> &initPrices,
                          int num_iterations) {

    vector<double> utility(num_bidders);

    for (int it = 0; it < num_iterations; ++it) {

        //in jeder iteration werden die preise des guts i auf die menge der preise,
        // die jeder bidder ausgegeben hat, gesetzt
        for (int j = 0; j < num_goods; ++j) {
            //TODO: verändertert, dass alle Budgets ausgegeben werden (oder alle items aufgebraucht)
            initPrices[j] = 0;
            for (int i = 0; i < bidders_PRD.size(); ++i)
                initPrices[j] += bidders_PRD[i].spent[j];

        }

        //update der valuations und spents pro bidder
        vector<vector<double>> update(bidders_PRD.size(), vector<double>(num_goods)); //
        for (int i = 0; i < bidders_PRD.size(); ++i) {
            for (int j = 0; j < num_goods; ++j) {
                update[i][j] = bidders_PRD[i].valuation[j] * bidders_PRD[i].spent[j] / initPrices[j];

            }
        }

        //new bid vector for next iteration
        for (int i = 0; i < bidders_PRD.size(); ++i) {
            for (int j = 0; j < num_goods; ++j) {
                bidders_PRD[i].spent[j] =  bidders_PRD[i].budget * update[i][j] / accumulate(update[i].begin(), update[i].end(), 0.0);

            }
        }



        if(it == (num_iterations - 1)){
            for (int j = 0; j < num_goods; ++j) {
                cout << initPrices[j] << "\n";

            }

        }



      /*//TODO: sind das korrekte Werte?
        if (it == (num_iterations - 1)) {
            for (int b = 0; b < num_bidders; ++b) {
                for (int i = 0; i < num_goods; ++i) {
                    utility[b] += double(bidders_PRD[b].valuation[i] * (bidders_PRD[b].spent[i] / initPrices[i]));
                }
                cout << "Utility (PR-D) für Bidder " << b << " ist: " << utility[b] << "\n";
                cout << "\n";
            }

        }*/


    }
    return initPrices;
}


/*
 * FIXME
 *  > Ergebnisse validieren !! Vor allem die Kosten !!!
 *  > Güter werden nicht komplett verkauft, d.h. keine market clearance
 *  > Utilities sind beim letzten Bidder am höchsten (?)
 *  > Budget wird (meist, außer bei sehr viele Bietern/Gütern) komplett aufgebraucht (!!)
 *     > initiale budget IST randomisiert
 *  > wie runden wir nun die fraktionalen Allocations ?
 *
 */


int main() {

    //generate #goods
    int num_goods;
    cout << "Number Goods: ";
    cin >> num_goods;

    //generate bidders
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

    /*//budget
    double budgetAgent;
    cout << "Budget of an agent: ";
    cin >> budgetAgent;*/

    // Anzahl Iterations pro 1 Handel (1x PR_Dynamics Algorithmus ausführen)
    int num_iterations;
    cout << "Number iterations: ";
    cin >> num_iterations;

    //quantity per item
    int quant;
    cout << "Quantity per item: ";
    cin >> quant;
    vector<double> quantItem(num_goods);
    for (int j = 0; j < num_goods; ++j) {
        quantItem[j] = quant;
    }

    //bidders vector for rounding
    vector<Bidder> bidders(num_bidders);

    vector<double> initBudget(num_bidders);

    double low_Budget = 1;
    double up_Budget = 10;

    double low_Val = 0;
    double up_Val = 11;

    for (int k = 0; k < num_bidders; ++k) {
        bidders[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        for (auto &v: bidders[k].valuation) v = (random_number(low_Val, up_Val));

        //budget
        initBudget[k] = (random_number(low_Budget, up_Budget));
        bidders[k].budget = initBudget[k];
        //= 1;
        //bidders[k].budget = budgetAgent;

        bidders[k].spent.resize(num_goods, bidders[0].budget / (double) num_goods);
    }


    //bidders vector for reference (PR-Dynamics algorithm)
    vector<Bidder> bidders_PRD(num_bidders);

    vector<double> initBudget_PRD(num_bidders);

    for (int k = 0; k < num_bidders; ++k) {
        bidders_PRD[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        for (auto &v: bidders_PRD[k].valuation) v = (random_number(low_Val, up_Val));

        //budget
        initBudget_PRD[k] = (random_number(low_Budget, up_Budget));
        bidders_PRD[k].budget = initBudget_PRD[k];
        //bidders[k].budget = budgetAgent;

        bidders_PRD[k].spent.resize(num_goods, bidders_PRD[0].budget / (double) num_goods);
    }




    //prices goods randomly initiated
    vector<double> initPrices(num_goods);

   /* for (int k = 0; k < num_goods; ++k) {
        //wichtig, dass hier ein spezifischer Bidder gewählt wird, da sonst Probleme bei num_goods > num_bidders und bidders[k].budget
        initPrices[k] = 1;
        //bidders[0].budget / num_goods;
    }*/

    //price update vector
    vector<vector<double>> update(bidders.size(), vector<double>(num_goods));


    vector<vector<double>> spendVec(num_bidders, vector<double>(num_goods));


    //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...
    vector<vector<myTuple>> mbbVec(num_bidders, vector<myTuple>(num_goods));

    for (int i = 0; i < num_bidders; ++i) {
        for (int k = 0; k < num_goods; ++k) {
            mbbVec[i][k].first = 1;
            mbbVec[i][k].second = k;
        }
    }


    //summe bisher spent per item
    vector<double> spendPerItem(num_goods);

    /*
    Funktionsaufrufe folgen hier:
    */

    //mbb graph
    vector<vector<myTuple>> SortedMbbVec = mbbGraph(num_bidders, num_goods, bidders, initPrices, mbbVec);


    //interests per good
    vector<int> interGood = interestsGood(num_bidders, num_goods, SortedMbbVec);


    //von Max utility und utility (im equilibrium sind diese gleich)
    vector<double> utility(num_bidders);
    vector<double> max_utility(num_bidders);


    /*
        Iterations folgen hier:
   */

    //rufe dann PR_D.cpp auf
    if(spendingRestriction == 0){

        //hier möchte ich funktion aufrufen aus PR_D.cpp
        PrDynamics_old(num_bidders, num_goods, bidders_PRD, initPrices,
                       num_iterations);

        cout << "\n";
        cout << "No restrictions. Therefore, we compute the PR-Dynamics Algorithm.";
        cout << "\n";
        exit(EXIT_FAILURE);



    }

    //für debugging
    //wiederholung für PR_D Algorithmus
    for (int it = 0; it < num_iterations; ++it) {


        ofstream myfile;
        myfile.open("data.txt", std::ios_base::app);

        //für budget printout (debugging)
        ofstream myfile2;
        myfile2.open("budget_final.txt", std::ios_base::app);



        //current price computation
        vector<double> newPrices = currentPrice(num_bidders, num_goods, bidders, initPrices, spendingRestriction,
                                                spendPerItem,
                                                quantItem, SortedMbbVec, interGood, spendVec, update, initBudget);


        //preisanpassung übernehmen für nächsten loop
        initPrices = newPrices;


        //falls alle Güter verkauft sind: FIXME: Market clearance
        if (accumulate(quantItem.begin(), quantItem.end(), 0.0) == 0.0) {
            cout << "Quantities are zero";
            cout << "\n";
            for (int i = 0; i < num_bidders; ++i) {
                if (bidders[i].budget < 0.01) {
                    bidders[i].budget = 0;
                }
                cout << "Budget bidder " << i << " is: " << bidders[i].budget << " (Budget zu " << (1-(bidders[i].budget/initBudget[i]))*100 << " % aufgebraucht!)";
                cout << "\n";
                myfile2 << "Budget bidder " << i << " is: " << bidders[i].budget << " (Budget zu " << (1-(bidders[i].budget/initBudget[i]))*100 << " % aufgebraucht!)";
                myfile2<< "\n";
            }
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < num_goods; ++i) {
            if (quantItem[i] < 0.01) {
                quantItem[i] = 0;
            }
        }


        //print for debugging
        cout << "\n";
        cout << "Iteration " << it << ":\n";
        for (int i = 0; i < num_goods; ++i) {
            cout << "Good " << i << " costs: " << newPrices[i] << "\n";
            //cout << "Good " << i << ": " << quantItem[i] << "\n";
        }
        cout << endl;



        //for debugging
        if (it == (num_iterations - 1)) {

            myfile << "Bidders: " << num_bidders << " Goods: " << num_goods << " Iterations: " << num_iterations
                   << " spending restriction: " << spendingRestriction << " quantity per item: " << quant << "\n";

          //Gleichgewichtspreise nach PR-Dynamics Algo
            initPrices = PrDynamics (num_bidders, num_goods, bidders_PRD, initPrices, num_iterations);


            cout << "\n";
            for (int i = 0; i < num_bidders; ++i) {
                if (bidders[i].budget < 0.01) {
                    bidders[i].budget = 0;
                }
            }

            //print for debugging
            cout << "\n";
            myfile << "\n";
            cout << "Final:\n";
            myfile << "Final:\n";
            for (int i = 0; i < num_goods; ++i) {
                cout << "Good " << i << " costs: " << newPrices[i] << " (PR-D: " << initPrices[i] << " )" << "\n";
                myfile << "Good " << i << " costs: " << newPrices[i] << " (PR-D: " << initPrices[i] << " )" << "\n";
            }
            cout << endl;
            myfile << endl;



            //print utility = valuation * menge (des guts) := bidders[i].valuation[j]*(spendVec[i][j]/newPrices[j])
            double utility = 0.0;

            //cout << "last, but not least";

            //print spending vector
            for (int i = 0; i < num_bidders; ++i) {
                cout << "Bidder " << i << " spends: " << "\n";
                myfile << "Bidder " << i << " spends: " << "\n";
                for (int j = 0; j < num_goods; ++j) {
                    cout << spendVec[i][j] << " | ";
                    myfile << spendVec[i][j] << " | ";
                }
                for (int j = 0; j < num_goods; ++j) {
                    utility += bidders[i].valuation[j] * (spendVec[i][j] / newPrices[j]);
                }
                cout << "\n";
                myfile << "\n";
                cout << "Utility: " << utility << "\n";
                myfile << "Utility: " << utility << "\n";
                cout << "Budget was: " << initBudget[i] << " (spend: " << 100*(1-(bidders[i].budget/initBudget[i])) << " %)" << "\n";
                myfile2 << "Budget was: " << initBudget[i] << " (spend: " << 100*(1-(bidders[i].budget/initBudget[i])) << " %)" << "\n";
                myfile << "Budget was: " << initBudget[i] << " (spend: " << 100*(1-(bidders[i].budget/initBudget[i])) << " %)"  << "\n";
                cout << "\n";
                myfile << "\n";
            }

            //wieviel bleibt pro Gut übrig:
            cout << "available items: \n";
            for (int j = 0; j < num_goods; ++j) {
                if (quantItem[j] < 0.01) {
                    quantItem[j] = 0;
                }
                cout << "Good " << j << " : " << quantItem[j] << "\n";
                myfile << "Good " << j << " : " << quantItem[j] << "\n";
            }

           myfile << "\n";


            /*  //for debuggin: testing if spendPerItem is met
             for (int j = 0; j < num_goods; ++j) {
                 cout << spendPerItem[j];
                 cout << "\n";
             }
 */

        }


    }


    return 0;
}

