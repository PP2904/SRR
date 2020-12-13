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
        for (auto &v: bidders[k].valuation) v = (random_number(4, 12));
        bidders[k].budget = 1;
    }

    //prices goods randomly initiated
    vector<double> prices(num_goods);

    for (int k = 0; k < num_goods; ++k) {
        //prices darf nicht 0 sein !!!
        prices[k] = double((random_number(6, 10))) / (random_number(1, 8));
        if (prices[k] == 0) {
            prices[k] = double((random_number(6, 10))) / (random_number(1, 8));
        }

    }

    //ich möchte einen bipartiten graphen erstellen (S. 6) und davon einen MBB graph und einen spending graph Q(x) ableiten
    // Dann möchte ich ein SR-equilibrium berechnen
    //danach möchte ich den SRR algo ausführen

    //berechnen MBB pro Bidder

    //alloc in mbb-graph
    vector<vector<double>> mbbAllocVec(num_bidders, vector<double>(num_goods));

    //tuple: (number good, mbb)
    typedef pair<double, int> mytuple;
    double mbb = 0.0;
    vector<mytuple> mbbVec(num_goods);

    for (int i = 0; i < num_bidders; ++i) {
        cout << "Bidder: " << i << "\n";
        for (int j = 0; j < num_goods; ++j) {
            mbb = bidders[i].valuation[j] / prices[j];
            mbbVec[j] = make_pair(mbb, j);
            // cout << "Bidder " << i << " MBB for Good " << j << " is: " << bidders[i].valuation[j]/prices[j] << "\n";
        }
        //sortiere mbb vector nach größe
        sort(mbbVec.begin(), mbbVec.end(), greater<>());

        //for debugging
        for (const mytuple &p: mbbVec) {
            cout << "Gut " << p.second << " hat MBB von: " << p.first << " und kostet " << prices[p.second] << "\n";
            mbbAllocVec[i][p.second] = p.first;
        }


    }


    //for debugging
    cout << "\n";
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            cout << mbbAllocVec[i][j] << " ";
        }
        cout << "\n";
    }

    //vector mit (größtem MBB,num_good) pro Bidder
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





    //logik for mbb-graph alloc


    //quanitity of item i intially = 1
    vector<double> quantItem(num_goods);
    for (int j = 0; j < num_goods; ++j) {
        quantItem[j] = 1;
    }

    cout << "\n";

    //alloc mbb items
    vector<vector<double>> mbbItemAllocVec(num_bidders, vector<double>(num_goods));

    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
            if (j == (greatestMBB[i].second)) {
                if (quantItem[j] != 0 && bidders[i].budget != 0) {
                    mbbItemAllocVec[i][j] = bidders[i].budget / prices[greatestMBB[i].second];
                    //item wurde vekauft und muss daher dezimiert werden
                    quantItem[j] = quantItem[j] - (bidders[i].budget / prices[greatestMBB[i].second]);
                    //stimm das mit dem budget abzug so?:                                                       
                    bidders[i].budget = bidders[i].budget - (bidders[i].budget / prices[greatestMBB[i].second]);
                    continue;
                }
                if (bidders[i].budget == 0) continue;

                if (quantItem[j] == 0) {

                    //TODO: gehe zu nächstem MBB-Item (leider hat greatestMBB nur die greatest MBBs...; müssen also mbbVec verwenden

                    mbbItemAllocVec[i][j] = 12.34;
                    continue;
                }

            } else {
                continue;
            }

        }
    }

    //debugging
    for (int i = 0; i < num_bidders; ++i) {
        for (int j = 0; j < num_goods; ++j) {
             cout <<  mbbItemAllocVec[i][j] << " ";
        }
        cout << "\n";
    }

    return 0;
}
