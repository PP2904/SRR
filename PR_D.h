//
// Created by Peter Pfrommer on 26.12.20.
//

#ifndef SRR_PR_D_H
#define SRR_PR_D_H

using namespace std;
#include <vector>

class Bidder {
public:
    vector<double> valuation; //was mir ein gut wert ist
    double budget;
    vector<double> spent; //für welches gut gibt bidder was aus (summe aller elem in spent ist budget)

    friend ostream &operator<<(ostream &os, const Bidder &b);
};



#endif //SRR_PR_D_H
