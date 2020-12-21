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
vector<vector<myTuple>> mbbGraph(int num_bidders, int num_goods, vector<Bidder> bidders, vector<double> newPrices) {


    double mbb;


    //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...
    vector<vector<myTuple>> mbbVec(num_bidders, vector<myTuple>(num_goods));

    //sorted vector for mbb
    vector<vector<myTuple>> SortedMbbVec(num_bidders, vector<myTuple>(num_goods));


    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {

            mbb = bidders[i].valuation[j] / newPrices[j];

            mbbVec[i][j] = make_pair(mbb, j);

        }

        //sortiere mbb vector nach größe bzgl. mbb wert (first value im tuple)


        sort(mbbVec[i].begin(), mbbVec[i].end(), greater<>());
        SortedMbbVec = mbbVec;
    }


    return SortedMbbVec;

}


// geht über den sortedVec und entscheidet pro Bidder, welche Güter er kauft;
// dabei darf maximal 1 GE auf jedes Gut verbraucht werden (in summe über alle Bidder)
vector<vector<double>> spendingGraph(int num_bidders, int num_goods, vector<Bidder> bidders, vector<double> newPrices,
                                     double spendingRestriction, vector<double> quantItem,
                                     vector<vector<myTuple>> SortedMbbVec) {

    //summe bisher spent per item
    vector<double> spendPerItem(num_goods);
    /*for (int j = 0; j < num_goods; ++j) {
        spendPerItem[j] = 0;
    }*/

    //spending vector  for spending graph Q(x);
    // Bidder xy: Gut1, Gut2, ...
    vector<vector<double>> spendVec(num_bidders, vector<double>(num_goods));

    /*for(int i = 0; i < num_bidders; ++i){
        for (int j = 0; j < num_goods; ++j) {
            spendVec[i][j] = 0.0;
        }
    }*/

    //spending graph erzeugen
    for (int iter = 0; iter < num_bidders; iter++) {
        //laufen durch den sortierten Vektor
        //spending-graph edges are a subset of the mbb-graph edges
        //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...
        for (const myTuple &p: SortedMbbVec[iter]) {

            //number des aktuellen Goods
            int numGood = p.second;

            //new share of good
            double newShare = (bidders[iter].budget / newPrices[p.second]);


            if (((spendPerItem[p.second] + (newShare * newPrices[p.second])) <= spendingRestriction) &&
                quantItem[p.second] != 0 && bidders[iter].budget != 0) {

                //spendVec wird erhöht durch neues share des Guts * price des guts
                spendVec[iter][p.second] = newShare * newPrices[p.second];

                //item wurde vekauft und muss daher dezimiert werden
                quantItem[p.second] -= newShare;

                //bisher für gut ausgegeben (overall agents)
                spendPerItem[p.second] += newShare * newPrices[p.second];


                // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                bidders[iter].budget -= newShare * newPrices[p.second];

                continue;
            }

            if (bidders[iter].budget == 0) {
                continue;
            }


            //if (((spendPerItem[p.second] + bidders[iter].budget / prices[p.second]) >= spendingRestriction)) { }


        }
    }


    return spendVec;
}

vector<int> interestsGood(int num_bidders, int num_goods, vector<vector<myTuple>> &SortedMbbVec) {

    //number of interested bidder in same good j (=interGood vector)
    vector<int> interGood(num_goods);


    //funktioniert das?
    for (int i = 0; i < num_bidders; ++i) {
        for (const myTuple &p: SortedMbbVec[i]) {

            if (p.first != 0.0) {
                interGood[p.second] += 1;
            }
        }

    }

    //for debugging
    cout << "\n";
    cout << "InterGood graph: \n";
    for (int j = 0; j < num_goods; ++j) {
            cout << interGood[j] << " ";
        }
        cout << "\n";



    return interGood;
}


vector<double> currentPrice(int num_bidders, int num_goods, vector<Bidder> bidders, vector<double> initPrices,
                            double spendingRestriction, vector<double> quantItem, vector<vector<myTuple>> &SortedMbbVec,
                            vector<int> &interGood) {
    //TODO: passen Preise schrittweise an

    vector<double> newPrices(num_goods);

    newPrices = initPrices;

    for (int k = 0; k < num_goods; ++k) {
        //wichtig, dass hier nur über goods gegangen wird, da sonst preise überschrieben werden
        //TODO: was ist wenn interGood = 0 ??
        // TODO: interGood richtig berechnet?

        if(interGood[k] > interGood[k]/2){
            newPrices[k] += newPrices[k] / interGood[k];
        }

        if(interGood[k] <= interGood[k]/2) {
            newPrices[k] -= newPrices[k] / interGood[k];
        }
    }


    //initPrices = newPrices;

    SortedMbbVec = mbbGraph(num_bidders, num_goods, bidders, newPrices);

    //TODO update spendVec bei jedem Loop
    vector<vector<double>> spendVec = spendingGraph(num_bidders, num_goods, bidders, newPrices, spendingRestriction, quantItem, SortedMbbVec);

    //for debugging
    cout << "\n";
    cout << "mbb-graph: \n";
    for (int i = 0; i < num_bidders; ++i) {
        cout << "Bidder " << i << ": " << "\n";
        //const iterator (läuft nur über das aktuelle mbbVec[i] = aktuelle Reihe i)
        for (const myTuple &p: SortedMbbVec[i]) {
            cout << "(" << p.second << "," << setprecision(3) << p.first << ")" << " ";
        }
        cout << "\n";
    }

    //for debugging
    cout << "\n";
    cout << "Spending graph: \n";
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout << spendVec[i][j] << " ";
        }
        cout << "\n";
    }






    return newPrices;

}


/*
 * Idee: 1) wir berechnen mbb graph, spending grap und die preise in separaten funktionen;
 *       2) jedes Mal, wenn wir den Preis verändern, geben wir mbb-graph und spending graph aus (für spending graph
 *              sind die Allokationen der Input := G(x) ) => Zusammenhang allocs und prices?
 *       3) der Preis für ein beliebtes Gut wird schrittweise erhöht, der für ein weniger beliebtes Gut erniedrigt (um wieviel?)
 *          => dadurch erhalten wir besser Infos über die Zuteilung von Güter zu Bietern
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

    //budget
    double budgetAgent;
    cout << "Budget of an agent: ";
    cin >> budgetAgent;

    //quantity item = 1
    vector<double> quantItem(num_goods);
    for (int j = 0; j < num_goods; ++j) {
        quantItem[j] = 1;
    }

    vector<Bidder> bidders(num_bidders);

    for (int k = 0; k < num_bidders; ++k) {
        bidders[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        for (int i = 0; i < num_goods; i++) {
            bidders[k].valuation[i] = (random_number(0, 6));
        }
        bidders[k].budget = budgetAgent;
    }


    //prices goods randomly initiated
    vector<double> initPrices(num_goods);

    for (int k = 0; k < num_goods; ++k) {
        //wichtig, dass hier ein spezifischer Bidder gewählt wird, da sonst Probleme bei num_goods > num_bidders und bidders[k].budget
        initPrices[k] = bidders[0].budget / num_goods;
    }


    //while( (accumulate(quantItem.begin(),quantItem.end(),0.0)) != 0.0 ) {
    /*
    Funktionsaufrufe folgen hier:
    */

    //mbb graph
    vector<vector<myTuple>> SortedMbbVec = mbbGraph(num_bidders, num_goods, bidders, initPrices);



    //für debugging
    for (int dur = 0; dur < 5; dur++) {

        //interests per good
        vector<int> interGood = interestsGood(num_bidders, num_goods, SortedMbbVec);


        //current price computation
        //TODO
        vector<double> newPrices = currentPrice(num_bidders, num_goods, bidders, initPrices, spendingRestriction,
                                                quantItem, SortedMbbVec, interGood);

        //preisanpassung übernehmen für nächsten loop
        initPrices = newPrices;




        //geht der funktionsaufruf?
        //vector<vector<double>> spendVec = spendingGraph(num_bidders, num_goods, bidders, newPrices, spendingRestriction, quantItem, SortedMbbVec);



    }


    return 0;
}

