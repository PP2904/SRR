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

//type def für: vector<vector<mytuple>> mbbVec
typedef pair<double, int> myTuple;


int random_number(int lb, int ub) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 engine(seed);
    return (engine() % (ub - lb + 1)) + lb;
}

//return mbbVec; vector<vector<myTuple>> mbbVec(num_bidders, vector<myTuple>(num_goods));
//bspw: bidder 1: (mbb, 0), (mbb,1), ...
vector<vector<myTuple>> mbbGraph(int num_bidders, int num_goods, vector<Bidder> &bidders, vector<double> &newPrices, vector<vector<myTuple>> &mbbVec) {


    //sorted vector for mbb
    vector<vector<myTuple>> SortedMbbVec(num_bidders, vector<myTuple>(num_goods));

    double mbb;




    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            //FIXME: macht die if-Abfrage Sinn?
           //if(mbbVec[i][j].first > 0) {

                mbb = bidders[i].valuation[j] / newPrices[j];
                if (mbb < 0.001) mbb = 0;

                mbbVec[i][j] = make_pair(mbb, j);
           //}
           continue;

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
                                     double spendingRestriction, vector<double> &quantItem,
                                     vector<vector<myTuple>> &SortedMbbVec, vector<vector<double>> spendVec) {

    //summe bisher spent per item
    vector<double> spendPerItem(num_goods);


    //spendVec = spending vector  for spending graph Q(x);
    // Bidder xy: Gut1, Gut2, ...


    //spending graph erzeugen
    for (int iter = 0; iter < num_bidders; iter++) {
        //laufen durch den sortierten Vektor
        //spending-graph edges are a subset of the mbb-graph edges
        //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...
        for (const myTuple &p: SortedMbbVec[iter]) {

            //number des aktuellen Goods
            int numGood = p.second;


            //TODO: ist hier noch ein Fehler?
            //new share of good
            double newShare = double(double(bidders[iter].budget/num_goods) / newPrices[p.second]);

            if (bidders[iter].budget == 0) {
                break;
            }


            if ((double(spendPerItem[p.second] + (newShare * newPrices[p.second])) <= spendingRestriction) &&
                    quantItem[p.second]-newShare >= 0.0 && bidders[iter].budget != 0.0) {


                //spendVec wird erhöht durch neues share des Guts * price des guts
                spendVec[iter][p.second] += newShare * newPrices[p.second];


                //item wurde vekauft und muss daher dezimiert werden
                quantItem[p.second] -= newShare;
                if(quantItem[p.second] < 0.001) quantItem[p.second] = 0.0;

                //bisher für gut ausgegeben (overall agents)
                spendPerItem[p.second] += newShare * newPrices[p.second];


             // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                bidders[iter].budget -= newShare * newPrices[p.second];

                break;
            }




        break;
        }
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
                            double spendingRestriction, vector<double> &quantItem, vector<vector<myTuple>> &SortedMbbVec,
                            vector<int> &interGood, vector<vector<double>> &spendVec, vector<vector<double>> &update) {

    ofstream myfile;
    myfile.open("data.txt", std::ios_base::app);

    //passen Preise schrittweise an
    vector<double> newPrices(num_goods);


    newPrices = initPrices;



      //PR-D preisanpassung; warum sinken die mbbs dennoch auf 0 mit iterations -> unendlich (in dem Fall schon 25 =) ??
      // Vermutung: die Preise werden immer geringer und dadurch das Verhältnis von valuation/preis ... )

    for (int k = 0; k < num_goods; ++k) {
        for (int i = 0; i < bidders.size(); ++i)
            newPrices[k] += bidders[i].spent[k];
    }

        //Problem ist, dass spent bereits im 2. Durchlauf = 0 ist. <- wird nicht upgedatet => update = 0 und dann teilen wir im
        //3. for loop hier durch update = 0 ...

        //erklärung: budget des bidders geht auf 0 (da er alles ausgegeben hat) und dann wird spent 0 => update = 0 und teilen durch 0 :(

        for (int i = 0; i < bidders.size(); ++i) {
            for (int j = 0; j < num_goods; ++j) {
                update[i][j] = bidders[i].valuation[j] * bidders[i].spent[j] / newPrices[j];

            }
        }

    for (int i = 0; i < bidders.size(); ++i) {
        for (int j = 0; j < num_goods; ++j) {
            bidders[i].spent[j] =
                    bidders[i].budget * update[i][j] / accumulate(update[i].begin(), update[i].end(), 0.0);

        }
    }


/* FUNKTIONSAUFRUFE */

    SortedMbbVec = mbbGraph(num_bidders, num_goods, bidders, newPrices, SortedMbbVec);

    interGood = interestsGood(num_bidders, num_goods, SortedMbbVec);

    //TODO update spendVec bei jedem Loop
    spendVec = spendingGraph(num_bidders, num_goods, bidders, newPrices, spendingRestriction, quantItem, SortedMbbVec, spendVec);





    //for debugging
    cout << "\n";
    cout << "mbb-graph: \n";
    for (int i = 0; i < num_bidders; ++i) {
        cout << "Bidder " << i << ": " << "\n";
        //const iterator (läuft nur über das aktuelle mbbVec[i] = aktuelle Reihe i)
        for (const myTuple &p: SortedMbbVec[i]) {
            if (i==0) {
                myfile << "(" << p.second << "," << setprecision(3) << p.first << ")" << " ";
            }
            cout << "(" << p.second << "," << setprecision(3) << p.first << ")" << " ";
        }
        cout << "\n";
        myfile << "\n";
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


/*
 * FIXME
 *  > Güter werden nicht komplett verkauft, d.h. keine market clearance
 *  > Utilities sind beim letzten Bidder am höchsten (?)
 *  > Budget wird komplett aufgebraucht (!!)
 *     > initiale budget IST randomisiert
 *  > wie runden wir nun die fraktionalen Allocations ?
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

    vector<Bidder> bidders(num_bidders);

    vector<double> initBudget(num_bidders);

    for (int k = 0; k < num_bidders; ++k) {
        bidders[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        for (auto &v: bidders[k].valuation) v = (random_number(0, 11));

        //budget
        initBudget[k] = (random_number(1, 7));
        bidders[k].budget = initBudget[k];
        //bidders[k].budget = budgetAgent;

        bidders[k].spent.resize(num_goods, bidders[0].budget / (double) num_goods);
    }


    //prices goods randomly initiated
    vector<double> initPrices(num_goods);

    for (int k = 0; k < num_goods; ++k) {
        //wichtig, dass hier ein spezifischer Bidder gewählt wird, da sonst Probleme bei num_goods > num_bidders und bidders[k].budget
        initPrices[k] = 1;
        //bidders[0].budget / num_goods;
    }

    //price update vector
    vector<vector<double>> update(bidders.size(), vector<double>(num_goods));


    vector<vector<double>> spendVec(num_bidders, vector<double>(num_goods));


    //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...
    vector<vector<myTuple>> mbbVec(num_bidders, vector<myTuple>(num_goods));

    for(int i = 0; i < num_bidders; ++i){
        for (int k = 0; k < num_goods; ++k) {
            mbbVec[i][k].first = 1;
            mbbVec[i][k].second = k;
        }
    }


    /*
    Funktionsaufrufe folgen hier:
    */

    //mbb graph
    vector<vector<myTuple>> SortedMbbVec = mbbGraph(num_bidders, num_goods, bidders, initPrices, mbbVec);


    //interests per good
    vector<int> interGood = interestsGood(num_bidders, num_goods, SortedMbbVec);

    //while( (accumulate(quantItem.begin(),quantItem.end(),0.0)) != 0.0 ) {

    //für debugging
    //wiederholung für PR_D Algorithmus
    for (int it = 0; it < num_iterations; ++it) {

        //current price computation
        vector<double> newPrices = currentPrice(num_bidders, num_goods, bidders, initPrices, spendingRestriction,
                                                quantItem, SortedMbbVec, interGood, spendVec, update);


        //preisanpassung übernehmen für nächsten loop
        initPrices = newPrices;


        if(accumulate(quantItem.begin(),quantItem.end(),0.0) == 0.0){
            cout << "Quantities are zero";
            exit(EXIT_FAILURE);
        }

        for(int i = 0; i < num_goods; ++i) {
            if (quantItem[i] < 0.01) {
                quantItem[i] = 0;
                continue;
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
       if(it == (num_iterations-1) ){

           //print utility
           double utility = 0.0;

           //cout << "last, but not least";

           //print spending vector
           for(int i = 0; i < num_bidders; ++i) {
               cout << "Bidder " << i << " : " << "\n";
               for (int j = 0; j < num_goods; ++j) {
                    cout << spendVec[i][j] << " | ";
               }
               for (int j = 0; j < num_goods; ++j) {
                   utility += newPrices[j]*spendVec[i][j];
               }
               cout << "\n";
               cout << "Utility: " << utility << "\n";
               cout << "Budget was: " << initBudget[i] << "\n";
               cout << "\n";
           }
           cout << endl;


        }


    }




    return 0;
}

