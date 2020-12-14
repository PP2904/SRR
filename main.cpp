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
        for(i=0; i < num_goods; i++){
               bidders[k].valuation[i] = (random_number(1, 12));
        }
        bidders[k].budget = 1;
    }

    //Idee 1: prices müssen zwischen 0 und 1 liegen, außer bei einem Gut

    //prices goods randomly initiated
    vector<double> prices(num_goods);

    for (int k = 0; k < num_goods; ++k) {

        //Test, ein Gut ist deutlich teurer als die andern
        if(k==(num_goods-1)){
          prices[k] = 1;
            continue;
        }

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

    //alloc in mbb-graph
    //vector<vector<double>> mbbAllocVec(num_bidders, vector<double>(num_goods));

    //tuple: (number good, mbb)
    typedef pair<double, int> mytuple;
    double mbb = 0.0;
    //vector Aufbau bspw: bidder 1: (mbb, 0), (mbb,1), ...
    vector<vector<mytuple>> mbbVec(num_bidders,vector<mytuple>(num_goods));

    for (int i = 0; i < num_bidders; ++i) {
        cout << "Bidder: " << i << "\n";
        for (int j = 0; j < num_goods; ++j) {

            //test: 1 bidder bewertet das künstlich veränderte teuerste Gut als besonders wertvoll
            if(i==(num_bidders-(num_bidders-1)) && j == (num_goods-1))
            {
                mbb = 100;
                mbbVec[i][j] = make_pair(mbb, j);
                continue;
            }
            
            mbb = bidders[i].valuation[j] / prices[j];
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
                cout << "(Element " << p.second << " mit MBB " << p.first  << " )" << " ";
            }
            //sehr wichtig!    
            break; 
        }
        cout << "\n";
    }

  /*  //vector mit (größtem MBB,num_good) pro Bidder
   vector<pair<double, int>> greatestMBB(num_bidders, make_pair(0, 0));

    cout << "\n";
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            if (mbbAllocVec[i][j] > greatestMBB[i].first) {
                greatestMBB[i].first = mbbAllocVec[i][j];
                greatestMBB[i].second = j;
            }
        }
        cout << "Für Bidder " << i << " greatest MBB " << "for good " << greatestMBB[i].second << " is: "
             << greatestMBB[i].first << "\n";

    }

    */

    cout << "\n";
    for (int i = 0; i < num_bidders; ++i) {
        cout << "Für Bidder " << i << " greatest MBB " << "\n";
        for (int j = 0; j < num_goods; ++j) {

            for (mytuple &p: mbbVec[i]) {

                cout  << " for good " << p.second << " is: " << p.first << "\n";
            }
            //sehr wichtig!
            break;

        }

    }



    //ACHTUNG:

    //logik for mbb-graph alloc:

      //every item is allocated fully
      //every agent spends all his budget
      //every agent spends budget only on MBB-items


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

    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {

            const mytuple &p = mbbVec[i][j];
            if (quantItem[j] != 0 && bidders[i].budget != 0) {
                //TODO: ist überhaupt noch genug für Bidder da?

                //allocation
                mbbItemAllocVec[i][j] = min(bidders[i].budget / prices[p.second], quantItem[j]);
                //spending
                spendVec[i][j] = mbbItemAllocVec[i][j] * prices[p.second];
                //item wurde vekauft und muss daher dezimiert werden
                quantItem[j] = quantItem[j] - mbbItemAllocVec[i][j];
                //stimm das mit dem budget abzug so?
                // (! //ACHTUNG: erst nach dem quantItem dezimiert ist, kann budget angepasst werde !)
                bidders[i].budget = bidders[i].budget - spendVec[i][j];
                continue;
            }
        }
    }



     // Printing the vectors

    //debugging - printing Alloc vector

    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
             cout <<  mbbItemAllocVec[i][j] << " ";
        }
        cout << "\n";
    }

    cout << "\n";

    //debugging - printing spending vector

    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout <<  spendVec[i][j] << " ";
        }
        cout << "\n";
    }


    //debugging - all goods allocated?
    cout << "\n";

    for (int j = 0; j < num_goods; ++j) {
        cout << quantItem[j] << " ";
    }



    return 0;
}
