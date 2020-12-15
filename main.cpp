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


struct Node {
    int data;
    struct Node *left;
    struct Node *right;

    // val is the key or the value that
    // has to be added to the data part
    Node(int val) {
        data = val;

        // Left and right child for node
        // will be initialized to null
        left = nullptr;
        right = nullptr;
    }
};

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


int main() {

    //generate #goods
    int num_goods;
    cout << "Number Goods: ";
    cin >> num_goods;


    //multiplier for valuation
    double i = 1.;

    //generate bidders with val, budget and spent_vec randomly
    int num_bidders;
    cout << "Number Bidders: ";
    cin >> num_bidders;


    //laut paper muss number bidders < num goods sein
    if (num_bidders > num_goods) {
        printf("Error number bidders larger than number goods");
        exit(EXIT_FAILURE);
    }


    vector<Bidder> bidders(num_bidders);

    for (int k = 0; k < num_bidders; ++k) {
        bidders[k].valuation.resize(num_goods);
        //valuation pro Gut und Bidder
        //for (auto &v: bidders[k].valuation){
        for (i = 0; i < num_goods; i++) {
            bidders[k].valuation[i] = (random_number(0, 8));
        }
        bidders[k].budget = 1;
    }

    //Idee 1: prices müssen zwischen 0 und 1 liegen, außer bei einem Gut

    //prices goods randomly initiated
    vector<double> prices(num_goods);

    for (int k = 0; k < num_goods; ++k) {

        /*//Test, ein Gut ist deutlich teurer als die andern
        if (k == (num_goods - 1)) {
            prices[k] = 1;
            continue;
        }*/

        prices[k] = double((random_number(6, 10))) / (random_number(10, 40));
        //prices darf nicht 0 sein !!!
        if (prices[k] == 0) {
            prices[k] = double((random_number(6, 10))) / (random_number(10, 50));
        }

    }

    //ich möchte einen bipartiten graphen erstellen (S. 6) und davon einen MBB graph und einen spending graph Q(x) ableiten
    // Dann möchte ich ein SR-equilibrium berechnen
    //danach möchte ich den SRR algo ausführen

    //berechnen MBB pro Bidder



    //number of interested bidder in same good j (=interGood vector)
    vector <int> interGood(num_goods);

    //initialize interGood vector
    for(int j = 0; j < num_goods; ++j){
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
            if(mbb != 0.0){
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

    //ACHTUNG:

    //logik for mbb-graph alloc:

    //every item is allocated fully
    //every agent spends all his budget
    //every agent spends budget only on MBB-items

    //maximum spend per item = 1
    vector<double> spendPerItem (num_goods);




    //quanitity of item i intially = 1
    vector<double> quantItem(num_goods);
    for (int j = 0; j < num_goods; ++j) {
        quantItem[j] = 1;
    }

    cout << "\n";

    //alloc mbb items
    vector<vector<double>> mbbItemAllocVec(num_bidders, vector<double>(num_goods));

    //spending vector  for spending graph Q(x)
    vector<vector<double>> spendVec(num_bidders, vector<double>(num_goods));

    //TODO: ich muss über mbbVec gehen und pro bidder das größte Gut zuweisen unter den Beschränkungen, die ich bereits verwendet habe
    // TODO: funktioniert ganz gut bisher, jedoch wieder die Schleife nur 1x ausgeführt; Ich benötige aber einen Durchlauf bis alle Items verkauft sind & jeder Bidder sein gesamtes Budget ausgegeben hat (!)


    while( (accumulate(quantItem.begin(),quantItem.end(),0)) != 0 ) {

        for (int iter = 0; iter < num_bidders; iter++) {

            for (const mytuple &p: mbbVec[iter]) {

                //number des aktuellen Goods
                int numGood = p.second;


                //Beschränkungen
                if (((spendPerItem[p.second] + min(bidders[iter].budget / prices[p.second], quantItem[p.second])) <=
                     1) && quantItem[p.second] != 0 && bidders[iter].budget != 0) {

                    //allocation; min <= ist überhaupt noch genug des Guts für den Bidder da?
                    //TODO: bidders[iter].budget / prices[p.second] = 1/0.28 ; quantItem[j] = 1;
                    //TODO: wir dürfen nur soviel zuweisen, wie auch einem Bidder zusteht; => wieviele bidder wollen das gut? => teile quantItem durch diese Anzahl und weise diese Menge maximal zu
                    mbbItemAllocVec[iter][p.second] = min(bidders[iter].budget / prices[p.second],
                                                          quantItem[p.second]); // /interGood[j]);
                    //spending
                    spendVec[iter][p.second] = mbbItemAllocVec[iter][p.second] * prices[p.second];
                    //spending per Item maximum of 1
                    spendPerItem[p.second] =
                            spendPerItem[p.second] + (mbbItemAllocVec[iter][p.second] * prices[p.second]);
                    //item wurde vekauft und muss daher dezimiert werden
                    quantItem[p.second] = quantItem[p.second] - mbbItemAllocVec[iter][p.second];
                    //stimm das mit dem budget abzug so?
                    // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                    bidders[iter].budget = bidders[iter].budget - spendVec[iter][p.second];

                    break;

                }

                if (((spendPerItem[p.second] + min(bidders[iter].budget / prices[p.second], quantItem[p.second])) >=
                     1)) {

                    //TODO springe zu nächste MBB des Bidders "iter"; also quasi nächste Tuple

                    const mytuple &p = mbbVec[iter][numGood];

                    //allocation; min <= ist überhaupt noch genug des Guts für den Bidder da?
                    //TODO: bidders[iter].budget / prices[p.second] = 1/0.28 ; quantItem[j] = 1;
                    //TODO: wir dürfen nur soviel zuweisen, wie auch einem Bidder zusteht; => wieviele bidder wollen das gut? => teile quantItem durch diese Anzahl und weise diese Menge maximal zu
                    mbbItemAllocVec[iter][p.second] = min(bidders[iter].budget / prices[p.second],
                                                          quantItem[p.second]); // /interGood[j]);
                    //spending
                    spendVec[iter][p.second] = mbbItemAllocVec[iter][p.second] * prices[p.second];
                    //spending per Item maximum of 1
                    spendPerItem[p.second] =
                            spendPerItem[p.second] + (mbbItemAllocVec[iter][p.second] * prices[p.second]);
                    //item wurde vekauft und muss daher dezimiert werden
                    quantItem[p.second] = quantItem[p.second] - mbbItemAllocVec[iter][p.second];
                    //stimm das mit dem budget abzug so?
                    // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                    bidders[iter].budget = bidders[iter].budget - spendVec[iter][p.second];

                    break;

                }


            }


        }

    }



/*    for (int j = 0; j < num_goods; ++j) {

        for (int iter = 0; iter < num_bidders; iter++) {

            const mytuple &p = mbbVec[iter][j];

            if (quantItem[j] != 0 && bidders[iter].budget != 0 && (spendPerItem[j] < 1)) {

                if((spendPerItem[j] + min(bidders[iter].budget / prices[p.second], quantItem[j]/interGood[j])) < 1 ){

                //allocation; min <= ist überhaupt noch genug des Guts für den Bidder da?
                //TODO: bidders[iter].budget / prices[p.second] = 1/0.28 ; quantItem[j] = 1;
                //TODO: wir dürfen nur soviel zuweisen, wie auch einem Bidder zusteht; => wieviele bidder wollen das gut? => teile quantItem durch diese Anzahl und weise diese Menge maximal zu
                mbbItemAllocVec[iter][j] = min(bidders[iter].budget / prices[p.second], quantItem[j]); // /interGood[j]);
                //spending
                spendVec[iter][j] = mbbItemAllocVec[iter][j] * prices[p.second];
                //spending per Item maximum of 1
                spendPerItem[j] = spendPerItem[j] + (mbbItemAllocVec[iter][j] * prices[p.second]);
                //item wurde vekauft und muss daher dezimiert werden
                quantItem[j] = quantItem[j] - mbbItemAllocVec[iter][j];
                //stimm das mit dem budget abzug so?
                // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                bidders[iter].budget = bidders[iter].budget - spendVec[iter][j];
                //continue;

                }
            }
        }
    }*/


    // Printing the vectors

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
        cout << quantItem[j] << " ";
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


    return 0;
}
