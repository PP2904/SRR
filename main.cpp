//
// Created by Peter Pfrommer on 19.12.20.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <iomanip>
#include <stdlib.h>
#include <utility>
#include <map>
#include <algorithm>
#include "mainPR_D.cpp"
#include "PR_D.h"
#include "rand_rounding.cpp"
#include "rand_rounding.h"

/* attention: hier wird der Code erklärt
 * Was passiert in SRR main:
 *   1) PR-Dynamics Preise werden pro iteration berechnet
 *   2) in jeder iteration wird der mbb-graph und der spending-graph berechnet
 *   3) diese graphen geben vor, welche Güter ein Bidder kaufen wird (möchte)
 *   4) hieraus ergeben sich also seine allocations
 *   5) diese allocations werden in rand_rounding.cpp gerundet und ausgegeben (mit der jeweiligen utility)
 *   6) falls restrictions = 0, dann wird der normale PR-Dynamics algo ausgeführt (das passiert in mainPR_D.cpp)
 *
 *
 *   Findings:
 *   - es werden im best case NUR 50% der items auf dem markt verkauft (available items aller Güter ist max. 50 % am Ende)
 *   - ich kann keine clearence des marktes (also auch kein Equilibrium erhalten), jedoch werden alle budgets
 *     voll ausgegeben
 *
 *
 */



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

//schon in PR_D.h vorhanden
/*ostream &operator<<(ostream &os, const Bidder &b) {
    for (int j = 0; j < b.spent.size(); ++j) {
        os << b.spent[j] << " ";
    }
    return os;
}*/

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
            //Attention: macht die if-Abfrage Sinn?
            //if(mbbVec[i][j].first > 0) {

            //valuations sind ints
            mbb = double(bidders[i].valuation[j]) / newPrices[j];
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
                                     vector<vector<myTuple>> &SortedMbbVec, vector<vector<double>> &spendVec,
                                     vector<vector<double>> &allocVec, vector<double> &allocVecOverall) {


    //spendVec = spending vector  for spending graph Q(x);
    // Bidder xy: Gut1, Gut2, ...



    //Attention: spendPerItem wird bei jeder iteration iter wieder auf 0 gesetzt, warum?
    // => weil keine Referenz (&) übergeben wurde von spendPerItem ...


    //das ist die neue Teil des Guts, welches der Bidder (welcher gerade an der Reihe ist) erwerben möchte
    double newShare = 0.0;

    //spending graph erzeugen
    for (int iter = 0; iter < num_bidders; iter++) {
        //laufen durch den sortierten Vektor
        //spending-graph edges are a subset of the mbb-graph edges
        //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...



        for (const myTuple &p: SortedMbbVec[iter]) {

            //number des aktuellen Goods
            int numGood = p.second;


            //ATTENTION: new share of good
            newShare = (bidders[iter].budget/double(num_goods))/newPrices[p.second];


            //ATTENTION: problem hier ist, dass bidder sein ganzes budget so ausgibt
            // -> sein budget ist 0 und dadurch ist spent = 0 -> update = 0 -> program exit ...
            //newShare = bidders[iter].budget/newPrices[p.second];




            if (bidders[iter].budget == 0) {
                break;
            }



            if ((double(spendPerItem[p.second] + (newShare * newPrices[p.second])) <= spendingRestriction) &&
                quantItem[p.second] - newShare >= 0.0 &&
                bidders[iter].budget != 0.0 &&
                //hier ist irgendwo noch ein Fehler, da die Bedingung nicht immer eingehalten wird ...
                //JA! DONE! wir haben keinen Vektor, der über alle Bidder summiert die allocation eines gutes überwacht!
                (allocVecOverall[p.second] + newShare) <= quantItem[p.second]) {

                allocVecOverall[p.second] += newShare;

                //spendVec wird erhöht durch neues share des Guts * price des guts
                spendVec[iter][p.second] += newShare * newPrices[p.second];

                allocVec[iter][p.second] += newShare;


                //item wurde vekauft und muss daher dezimiert werden
                quantItem[p.second] -= newShare;
                if (quantItem[p.second] < 0.1) quantItem[p.second] = 0.0;

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


vector<double> currentPrice(int num_bidders, int num_goods, vector<Bidder> &bidders, vector<double> &initPrices,
                            double spendingRestriction, vector<double> &spendPerItem, vector<double> &quantItem,
                            vector<vector<myTuple>> &SortedMbbVec,
                            vector<int> &interGood,
                            vector<vector<double>> &spendVec, vector<vector<double>> &update,
                            vector<double> &initBudget, vector<vector<double>> &allocVec,
                            int quant, int num_iterations, vector<double> &allocVecOverall) {

    ofstream myfile;
    myfile.open("data.txt", std::ios_base::app);

    ofstream rdResult;
    rdResult.open("roundedResult.txt", std::ios_base::app);

    //passen Preise schrittweise an
    vector<double> newPrices(num_goods);

    vector<double> MaxUtility(num_bidders);
    double utility = 0.0;


    newPrices = initPrices;



    //PR-D preisanpassung; warum sinken die mbbs dennoch auf 0 mit iterations -> unendlich (in dem Fall schon 25 =) ??
    // Vermutung: die Preise werden immer geringer und dadurch das Verhältnis von valuation/preis ... )

    for (int k = 0; k < num_goods; ++k) {
        //Attention: verändert, dass alle Budgets ausgegeben werden (oder alle items aufgebraucht) - entweder oder (!)
        // Attention: Aussage von oben stimmt nicht mehr (6.1.21)
        newPrices[k] = 0;
        for (int i = 0; i < bidders.size(); ++i) {
            newPrices[k] += bidders[i].spent[k];
        }
    }



    //Problem ist, dass spent bereits im 2. Durchlauf = 0 ist. <- wird nicht upgedatet => update = 0 und dann teilen wir im
    //3. for loop hier durch update = 0 ...

    //erklärung: budget des bidders geht auf 0 (da er alles ausgegeben hat) und dann wird spent 0 => update = 0 und teilen durch 0 :(

    for (int i = 0; i < bidders.size(); ++i) {
        for (int j = 0; j < num_goods; ++j) {
            update[i][j] = double(bidders[i].valuation[j] * bidders[i].spent[j]) / newPrices[j];

            //Attention for debugging
            if(isnan(update[i][j])){
                cout << "update isnan";
                rdResult << "update isnan";
                //Trick siebzehn
                newPrices[1] = 999;
                return newPrices;
                //exit(EXIT_FAILURE);
            }
        }
    }


    for (int i = 0; i < bidders.size(); ++i) {
        for (int j = 0; j < num_goods; ++j) {

            ofstream rdResult;
            rdResult.open("roundedResult.txt", std::ios_base::app);

            //Attention: Problem: irgendwann ist der update vector 0, weil spent-vector sehr klein ist und daher
            // nur noch wenig hinzukommt;


            //Attention: hier wird eingestiegen, wenn update vector 0 (fast 0) ist (!)
            if (accumulate(update[i].begin(), update[i].end(), 0.0) <= 0.001) {


                //debugging
                /*cout << "update zero";
                exit(EXIT_FAILURE);*/


                cout << endl;
                //wieviel bleibt pro Gut übrig:
                cout << "available items: \n";
                for (int j = 0; j < num_goods; ++j) {
                    if (quantItem[j] < 0.1) {
                        quantItem[j] = 0;
                    }
                    cout << "Good " << j << " : " << quantItem[j] << "\n";
                    myfile << "Good " << j << " : " << quantItem[j] << "\n";
                }
                for (int i = 0; i < bidders.size(); ++i) {
                    cout << "(Bidder " << i << " (spend: " << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)"
                         << "\n";
                    //rdResult << "Bidder " << i << " Budget was: " << initBudget[i] << " (spend: "
                    //  << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";

                    //Attention: hier printen wir aktuelles budget und die Ausgaben in % an
                    rdResult << initBudget[i] << " | " << 100 * (1 - (bidders[i].budget / initBudget[i])) << "\n";

                    for (int j = 0; j < num_goods; ++j) {
                        utility += bidders[i].valuation[j] * allocVec[i][j];
                        if (j == (num_goods - 1)) {
                            MaxUtility[i] = utility;
                        }
                    }


                }

                rdResult << "\n";

                //for debugging
                for (int j = 0; j < num_goods; ++j) {
                    cout << "\n";
                    cout << "Für Good " << j << " wurden " << spendPerItem[j] << " Geldeinheiten ausgegeben" << "\n";
                }
                cout << "\n";


                //Attention: warum hier aufruf von roundingSRE ??
                //print to file "myfile" findet in rand_rounding.cpp statt
                //Spending Restricted Equilibrium (SRE)
                roundingSRE(num_bidders, num_goods, allocVec, quant, bidders, MaxUtility, spendingRestriction,
                            num_iterations, newPrices);

                cout << "\n";
                cout << "update vector is zero";
                rdResult << "\n";
                rdResult << "update vector is zero";
                rdResult << "\n";
                //Trick siebzehn
                newPrices[1] = 999;
                return newPrices;
                //exit(EXIT_FAILURE);
                //break;
            }
            bidders[i].spent[j] =
                    (bidders[i].budget * update[i][j]) / accumulate(update[i].begin(), update[i].end(), 0.0);

        }
    }





/* FUNKTIONSAUFRUFE */



    SortedMbbVec = mbbGraph(num_bidders, num_goods, bidders, newPrices, SortedMbbVec);

    interGood = interestsGood(num_bidders, num_goods, SortedMbbVec);


    spendVec = spendingGraph(num_bidders, num_goods, bidders, newPrices, spendingRestriction, spendPerItem, quantItem,
                             SortedMbbVec,
                             spendVec, allocVec, allocVecOverall);


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
            //Attention: verändert, dass alle Budgets ausgegeben werden (oder alle items aufgebraucht)
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
                bidders_PRD[i].spent[j] =
                        bidders_PRD[i].budget * update[i][j] / accumulate(update[i].begin(), update[i].end(), 0.0);

            }
        }


        //debugging
        /*if(it == (num_iterations - 1)){
            for (int j = 0; j < num_goods; ++j) {
                cout << initPrices[j] << "\n";

            }

        }*/


    }
    return initPrices;
}


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

    //Attention: iterations overall (Gesamtexperiment wdh)
    int num_iter_exp;
    cout << "Number iterations overall: ";
    cin >> num_iter_exp;

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

    double low_Budget = 10;
    double up_Budget = 30;

    double low_Val = 0;
    double up_Val = 15;

    for (int k = 0; k < num_bidders; ++k) {
        bidders[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        for (auto &v: bidders[k].valuation) v = (random_number(low_Val, up_Val));

        //budget
        initBudget[k] = (random_number(low_Budget, up_Budget));
        bidders[k].budget = initBudget[k];
        //= 1;
        //bidders[k].budget = budgetAgent;

        //attention 8.1.21: *10
        bidders[k].spent.resize(num_goods, bidders[0].budget / (double) num_goods*10);
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

        //attention 8.1.21: *10
        bidders_PRD[k].spent.resize(num_goods, bidders_PRD[0].budget / (double) num_goods*10);
    }




    //prices of goods randomly initiated
    vector<double> initPrices(num_goods);


    //price update vector
    vector<vector<double>> update(bidders.size(), vector<double>(num_goods));


    vector<vector<double>> spendVec(num_bidders, vector<double>(num_goods));

    //alloc vector for one bidder
    vector<vector<double>> allocVec(num_bidders, vector<double>(num_goods));

    // alloc vector overall bidders
    vector<double> allocVecOverall(num_goods, 0.0);


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
    if (spendingRestriction == 0) {

        //hier möchte ich funktion aufrufen aus PR_D.cpp
        PrDynamicsRestZero(num_bidders, num_goods, bidders_PRD, initPrices,
                           num_iterations);

        cout << "\n";
        cout << "No restrictions. Therefore, we compute the PR-Dynamics Algorithm.";
        cout << "\n";
        exit(EXIT_FAILURE);


    }


    /* ATTENTION
     * ATTENTION: Hier beginnen die Iterations (außen: für Gesamtexperiment, innen: für PR-Dynamics)
     * ATTENTION
     */


    //FOR SCHLEIFE FÜR ANZAHL WIEDERHOLUNGEN DES GESAMTEXPERIMENTS
    for (int iter = 0; iter < num_iter_exp; iter++) {

        //Attention: zurücksetzen der Werte pro Gesamtexperiment-Wdh.
        if (iter > 0) {

            //Attention: gleiche Eigenschaften wie zu beginn wiederherstellen
            for (int k = 0; k < num_bidders; ++k) {
                //bidders[k].valuation.resize(num_goods);
                //valuation pro Gut und Bidder
                //for (auto &v: bidders[k].valuation) v = (random_number(low_Val, up_Val));

                //budget
                //initBudget[k] = (random_number(low_Budget, up_Budget));
                bidders[k].budget = initBudget[k];


                /*//attention 8.1.21: *10
                bidders[k].spent.resize(num_goods, bidders[0].budget / (double) num_goods*10);
*/

                //attention 8.1.21:
                //setzen spendVec pro Overall Iteration wieder auf 0
                for (int j = 0; j < num_goods; ++j){
                    spendVec[k][j] = 0.0;
                    allocVec[k][j] = 0.0;
                    update [k][j] = 0.0;

                    bidders[k].spent[j] = bidders[0].budget / (double) num_goods*10;



                }


                for (int j = 0; j < num_goods; ++j){
                SortedMbbVec[k][j].first = 0.0;
                SortedMbbVec[k][j].second = 0;
                }


            }

            //zurücksetzen
            for (int j = 0; j < num_goods; ++j){
                spendPerItem[j] = 0.0;
                quantItem[j] = quant;
            }

            //prices wieder zurücksetzen für currentPrices()
            for (int j = 0; j < num_goods; ++j){
                initPrices[j]=0.0;
                allocVecOverall[j] = 0.0;
                interGood[j] = 0.0;
            }



            /*
             * HIER KÖNNEN MEHRERE EIGENSCHAFTEN STEHEN
             */

/*
            //spending restriction
            cout << "What is the spending restriction?: ";
            cin >> spendingRestriction;*/




            //quantity per item
            cout << "Quantity per item: ";
            cin >> quant;
            for (int j = 0; j < num_goods; ++j) {
                quantItem[j] = quant;
            }

        }


        //für debugging
        //wiederholung für PR_D Algorithmus
        for (int it = 0; it < num_iterations; ++it) {


            ofstream myfile;
            myfile.open("data.txt", std::ios_base::app);

            //für budget printout (debugging) FALLS quantities nicht 0 (!)
            ofstream myfile2;
            myfile2.open("resultNotZeroQuant.txt", std::ios_base::app);

            //für budget printout (debugging)
            ofstream myfileZEROQuant;
            myfileZEROQuant.open("resultsZeroQuant.txt", std::ios_base::app);

            ofstream rdResult;
            rdResult.open("roundedResult.txt", std::ios_base::app);



            //current price computation
            vector<double> newPrices = currentPrice(num_bidders, num_goods, bidders, initPrices, spendingRestriction,
                                                    spendPerItem,
                                                    quantItem, SortedMbbVec, interGood, spendVec, update, initBudget,
                                                    allocVec, quant, num_iterations, allocVecOverall);

            //falls update zero in currentPrice funktion
            if(newPrices[1] == 999) break;

            //preisanpassung übernehmen für nächsten loop
            initPrices = newPrices;


            /*
             *
             * QUANTITIES are ZERO
             *
             * */


            //Attention: falls alle Güter verkauft sind: Market clearance
            if (accumulate(quantItem.begin(), quantItem.end(), 0.0) == 0.0) {


                vector<double> utilityRound(num_bidders);
                double utilityQuantZero = 0.0;


                myfileZEROQuant << "Quantities are zero";
                cout << "\n";
                myfileZEROQuant << "\n";


                for (int i = 0; i < num_bidders; ++i) {
                    //davor war < 0.01
                    if (bidders[i].budget < 0.1) {
                        bidders[i].budget = 0;
                    }
                    cout << "Budget bidder " << i << " is: " << bidders[i].budget << " (Budget zu "
                         << (1 - (bidders[i].budget / initBudget[i])) * 100 << " % aufgebraucht!)";
                    cout << "\n";
                    myfileZEROQuant << "Budget bidder " << i << " is: " << bidders[i].budget << " (Budget zu "
                                    << (1 - (bidders[i].budget / initBudget[i])) * 100 << " % aufgebraucht!)";
                    myfileZEROQuant << "\n";
                    cout << "Bidder " << i << " spends: \n";
                    for (int j = 0; j < num_goods; ++j) {
                        cout << spendVec[i][j] << " | ";
                    }
                    cout << "\n";
                    cout << "Bidder " << i << " allocation: \n";
                    for (int j = 0; j < num_goods; ++j) {

                        //alloc of goods for each bidder (findet in spendingGraph statt)
                        //allocVec[i][j] = double(spendVec[i][j] / newPrices[j]);

                        cout << allocVec[i][j] << " | ";
                        myfileZEROQuant << allocVec[i][j] << " | ";
                    }



                    /* }


                     for (int i = 0; i < num_bidders; ++i) {*/
                    for (int j = 0; j < num_goods; ++j) {
                        //utilityQuantZero += bidders[i].valuation[j] * (spendVec[i][j] / newPrices[j]);
                        utilityQuantZero += bidders[i].valuation[j] * allocVec[i][j];
                        if (j == (num_goods - 1)) {
                            utilityRound[i] = utilityQuantZero;
                        }
                    }

                    cout << "\n";
                    myfile << "\n";
                    cout << "Utility: " << utilityQuantZero << "\n";
                    myfile << "Utility: " << utilityQuantZero << "\n";
                    myfileZEROQuant << "Utility: " << utilityQuantZero << "\n";
                    cout << "Budget was: " << initBudget[i] << " (spend: "
                         << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";
                    myfileZEROQuant << "Budget was: " << initBudget[i] << " (spend: "
                                    << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";
                    myfile << "Budget was: " << initBudget[i] << " (spend: "
                           << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";
                    /*rdResult << "Bidder " << i << " Budget was: " << initBudget[i] << " (spend: "
                             << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";*/

                    //Attention: hier printen wir aktuelles budget und die Ausgaben in % aus
                    rdResult << initBudget[i] << " | " << 100 * (1 - (bidders[i].budget / initBudget[i])) << "\n";


                    cout << "\n";
                    myfile << "\n";
                    myfileZEROQuant << "\n";

                }
                rdResult << "\n";

                //falls quantities zero
                int flag;

                roundingSRE(num_bidders, num_goods, allocVec, quant, bidders, utilityRound, spendingRestriction,
                            num_iterations, newPrices);

                cout << "\n";
                cout << "Quantities are zero";
                rdResult << "\n";
                rdResult << "Quantities are zero";
                rdResult << "\n";
                //exit(EXIT_FAILURE);
                break;
            }

            for (int i = 0; i < num_goods; ++i) {
                if (quantItem[i] < 0.1) {
                    quantItem[i] = 0;
                }
            }


            //print for debugging
            cout << "\n";
            cout << "Iteration " << it << ":\n";
            for (int i = 0; i < num_goods; ++i) {
                cout << "Good " << i << " costs: " << newPrices[i] << "\n";
                cout << "Quantity Good " << i << ": " << quantItem[i] << "\n";
            }
            cout << endl;



            /*
            *
            * QUANTITIES are NOT ZERO
            *
            * */


            //Attention: hier steigen wir ein, wenn alles ohne Probleme durchläuft (update vector bleibt >0 und quantities werden nicht 0)
            //for debugging
            if (it == (num_iterations - 1)) {

                myfile << "Bidders: " << num_bidders << " Goods: " << num_goods << " Iterations: " << num_iterations
                       << " spending restriction: " << spendingRestriction << " quantity per item: " << quant << "\n";

                //Gleichgewichtspreise nach PR-Dynamics Algo
                initPrices = PrDynamics(num_bidders, num_goods, bidders_PRD, initPrices, num_iterations);


                cout << "\n";
                for (int i = 0; i < num_bidders; ++i) {
                    if (bidders[i].budget < 0.1) {
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
                    rdResult << newPrices[i] << " | " << initPrices[i] << "\n";
                }
                cout << endl;
                myfile << endl;



                //print utility = valuation * menge (des guts) := bidders[i].valuation[j]*(spendVec[i][j]/newPrices[j])
                double utility = 0.0;

                //vector der utility pro Bidder für rand_rounding
                // nicht das gleiche wie max_utility
                vector<double> MaxUtility(num_bidders);

                //cout << "last, but not least";

                //print spending vector
                //attention: Bidder 0 spends 32.5 | 8.89 | 44.9 | 5.21 | 11.1 | 13.5 |, but budget is: 31 only (?!)
                // gelöst 8.1.21
                for (int i = 0; i < num_bidders; ++i) {
                    utility = 0.0;
                    cout << "Bidder " << i << " spends: " << "\n";
                    myfile << "Bidder " << i << " spends: " << "\n";
                    for (int j = 0; j < num_goods; ++j) {
                        cout << spendVec[i][j] << " | ";
                        myfile << spendVec[i][j] << " | ";
                    }
                    cout << "\n";
                    myfile << "\n";
                    cout << "Bidder " << i << " allocation: " << "\n";
                    myfile << "Bidder " << i << " allocation: " << "\n";
                    for (int j = 0; j < num_goods; ++j) {


                        //alloc of goods for each bidder (findet in spendingGraph statt)
                        //allocVec[i][j] = double(spendVec[i][j] / newPrices[j]);

                        cout << allocVec[i][j] << " | ";
                        myfile << allocVec[i][j] << " | ";

                    }
                    /* }


                     for (int i = 0; i < num_bidders; ++i) {*/
                    for (int j = 0; j < num_goods; ++j) {
                        //utility += bidders[i].valuation[j] * (spendVec[i][j] / newPrices[j]);
                        utility += bidders[i].valuation[j] * allocVec[i][j];
                        if (j == (num_goods - 1)) {
                            MaxUtility[i] = utility;
                        }
                    }


                    cout << "\n";
                    myfile << "\n";
                    cout << "Utility: " << utility << "\n";
                    myfile << "Utility: " << utility << "\n";
                    cout << "Budget was: " << initBudget[i] << " (spend: "
                         << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";
                    myfile2 << "Budget was: " << initBudget[i] << " (spend: "
                            << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";
                    myfile << "Budget was: " << initBudget[i] << " (spend: "
                           << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";
                    /*rdResult << "Bidder " << i << " Budget was: " << initBudget[i] << " (spend: "
                             << 100 * (1 - (bidders[i].budget / initBudget[i])) << " %)" << "\n";*/
                    //Attention: hier schreiben wir in die Datei
                    rdResult << initBudget[i] << " | " << 100 * (1 - (bidders[i].budget / initBudget[i])) << "\n";


                    cout << "\n";
                    myfile << "\n";
                }

                rdResult << "\n";

                //wieviel bleibt pro Gut übrig:
                cout << "available items: \n";
                rdResult << "available items: \n";
                for (int j = 0; j < num_goods; ++j) {
                    if (quantItem[j] < 0.1) {
                        quantItem[j] = 0;
                    }
                    cout << "Good " << j << " : " << quantItem[j] << "\n";
                    myfile << "Good " << j << " : " << quantItem[j] << "\n";
                    rdResult << quantItem[j] << " | ";
                }

                rdResult << "\n";
                rdResult << "iter: " << "\n";
                rdResult << iter << "\n";
                rdResult << "\n";


                myfile << "\n";
                cout << "\n";



                //MaxUtility ist die utility pro bidder (in rand iterative rounding ist das max_utility

                //print to file "myfile" findet in rand_rounding.cpp statt
                //Spending Restricted Equilibrium (SRE)
                roundingSRE(num_bidders, num_goods, allocVec, quant, bidders, MaxUtility, spendingRestriction,
                            num_iterations, newPrices);


            }


        }

        //Attention: für overall iterations
    }


    return 0;


}